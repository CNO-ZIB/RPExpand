function out = f(sc, varargin)
%F Solves scattering problems at complex frequencies
%   tags = F(sc, contours) This function is called when an instance sc of
%   class Scattering is called with (). The input argument contours is 
%   expected to be a cell array of size=(1,:) containing vectors of
%   complex frequencies. It solves the scattering problems along the given
%   contours. If the input is an empty cell {} is returned. Otherwise a
%   cell array of the same shape as contours containing tags corresponding
%   to solutions of the scattering problems.
%
%   v = F(sc, tags, quantity) The specified quantity is evaluated in a post
%   process from solutions of the linear systems along the contours, which 
%   are represented by tags = f(sc, contours). Note, that for quadratic
%   quantities solutions for two different frequencies are required.
%
%   See also SCATTERING, RIESZPROJECTION

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

if isempty(varargin{1}), out = {}; return; end
if length(varargin)==1
    out = solveScatteringProblems(sc,varargin{1});
else
    out = sc.postprocess(varargin{:});
end
end

function tags = solveScatteringProblems(sc,contours)
% check if a resource has been set up for computation
if isempty(jcmwave_daemon_resource_info)
    error('Please start a daemon');
end

% select or create resultbag
if isempty(contours), tags=0; return; end
isRef = length(contours)==1 && isreal(contours{1});
if isRef, scDir = sc.referenceDir; else, scDir = sc.scatteringDir; end
index = isRef+1;
if isempty(sc.resultbags{index})
    bag = sc.bags{index};
    filename = [sc.resultbagDir filesep bag];
    relevantfields = fieldnames(sc.keys).';
    relevantfields(strcmp('relPerm',relevantfields)) = [];
    resultbag = jcmwave_resultbag(filename,relevantfields);
    sc.resultbags{index} = resultbag;
else
    resultbag = sc.resultbags{index};
end

keys = sc.keys;

nContours = length(contours);
nPoints = sum(cellfun(@length, contours));

if ~nPoints, tags = cell(1,nContours); return; end

jcmpSc = sc.projectFile;

% check if field exports are required
pointList = [];
pev = {'DipoleEmission' 'PointEvaluation'};
pev = pev(cellfun(@(x)sc.quantities.(x).logical(3),pev));
if ~isempty(pev) && sc.quantities.(pev{1}).logical(3)
    pointList = sc.quantities.(pev{1}).position;
end

% check if the postprocess has been removed and if derivatives are required
id = fopen(jcmpSc, 'r');
jcmpt = fread(id,'*char')';
fclose(id);
dx = regexp(jcmpt,'(?<=Derivatives\s*{.+Order\s*=\s*)\S+','match','once');
if ~isempty(dx) 
    if dx(1)=='%'
        dx = keys.(regexp(dx,'(?<=\()\w+','match','once')); 
    else
        dx = dx(1)=='1';
    end
else
    dx = false;
end

sc.derivatives = dx;

s = regexp(jcmpt,sc.exportPoints(1:29),'once');
if ~isempty(s), jcmpt(s:end) = []; end
if isempty(pointList) && ~isempty(s)
    id = fopen(jcmpSc,'w');
    fwrite(id,jcmpt,'char');
    fclose(id);
% append a post process to the project file if needed
elseif ~isempty(pointList)
    if isRef && strcmp(pev{1},'PointEvaluation')
        sc.quantities.(pev{1}).reference.resultbag = resultbag;
    elseif ~isRef
        sc.quantities.(pev{1}).resultbag = resultbag;
    end
    exportPoints = sc.exportPoints;
    if polyval(sc.version>[5;0;6],2)>4 % 6.0.7 or newer
        keys.localFieldOnly = sprintf('LocalResponseField = %s','yes');
    else
        keys.localFieldOnly = sprintf('AddSingularFields = %s','no');
    end
    id = fopen(jcmpSc,'w');
    fwrite(id,[jcmpt sprintf(exportPoints)],'char');
    fclose(id);
    resdir = [sc.projectFile(length(sc.projectDir)+2:end-6) '_results'];
    keys.fPath = [resdir filesep 'fieldbag.jcm'];
    keys.oPath = [resdir filesep 'points.jcm'];
    keys.px = pointList(:,1);
    keys.py = pointList(:,2);
    keys.pz = pointList(:,3);
    keys.quantity = 'ElectricFieldStrength';
    keys.format = 'JCM';
    keys.keepDerivatives = 'no';
    if dx, keys.format = 'JCM-ASCII'; keys.keepDerivatives = 'yes'; end
end

jcmpSc = jcmpSc(1:end-1);
job_ids = zeros(1,nPoints,'int32');

% evaluate first point to get common PMLs if there is no existing log file
logFile = sc.pml; idx = 0; waitForLogfile = false;
if isempty(logFile)
    logFile = [sc.workingDir filesep 'pml.log']; waitForLogfile = true;
end
if ~exist(logFile, 'file') && ~isempty(contours{1})
    keys.pml = sprintf('LogFile = "%s"', logFile);
    keys.omega = complex(contours{1}(1));
    if isfield(keys, 'relPerm') && isa(keys.relPerm, 'function_handle')
        keys.relPermittivity = keys.relPerm(keys.omega);
    end
    wdir = [scDir filesep resultbag.get_tag(keys)];
    job_ids(1) = jcmwave_solve(jcmpSc,keys,resultbag,'workingdir',wdir);
    if ~job_ids(1), waitForLogfile = false; end
    idx = idx+1;
end

% evaluate the remaining points
t0 = tic;
while ~exist(logFile, 'file') && waitForLogfile
    if toc(t0)>60
        jcmwave_daemon_wait(job_ids(1), resultbag);
        logs = resultbag.get_log(keys);
        fprintf(1,logs.Log.Out);
        error(['\nThe file ''pml.log'' is missing. Probably the '...
               'PML section in your jcmpt file is missing']);
    end
    pause(1); 
end
tags = cell(1,nContours); keys.pml = sprintf('InputFile = "%s"', logFile); 
for it1 = 1:nContours
    tags_it1 = cell(1,numel(contours{it1}));
    for it2 = 1:numel(contours{it1})
        keys.omega = complex(contours{it1}(it2));
        if isfield(keys, 'relPerm')&&isa(keys.relPerm, 'function_handle')
            keys.relPermittivity = keys.relPerm(keys.omega);
        end
        tags_it1{it2} = resultbag.get_tag(keys);
        wdir = [scDir filesep tags_it1{it2}]; idx = idx+1;
        job_ids(idx) = jcmwave_solve(jcmpSc,keys,resultbag, ...
            'workingdir',wdir);
    end
    tags{it1} = tags_it1;
end

% delete temporary jcmpt file
if ~isempty(pointList)
    id = fopen([jcmpSc 't'],'w');
    fwrite(id,jcmpt,'char');
    fclose(id);
end
jcmwave_daemon_wait(job_ids, resultbag);
end
