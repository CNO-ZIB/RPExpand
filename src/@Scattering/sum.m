function out = sum(sc,tags,weights,omega,outName)
%SUM Compute a weighted sum of field bags or field values
%   out = SUM(sc, tags, weights, omega0, outName) returns the weighted sum 
%   of the fields represented by the tags. The post process
%   <a href="matlab:web(...
%   ['https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/' ...
%   '721244647add8d817b99d2e4661c6ca8.html'], '-browser')"
%   >Superposition</a> of JCMsuite is used. If specified the outNames
%   determine the output filename. Omega is the frequency of the result.
%   The job ids and the corresponding tags are returned.
%
%   see also RieszProjection, jcmwave_resultbag,
%            jcmwave_solve, jcmwave_daemon_wait

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

a = filesep; pDir = [sc.workingDir a 'superpositions'];

% check if a resource has been set up for computation
if isempty(jcmwave_daemon_resource_info)
    error('Please start a daemon');
end

if nargin<5, outName = ''; end

% get resultbag with scattering resultss
if exist(tags{1},'file') || isempty(sc.resultbags{2}) || ...
        ~isfield(sc.resultbags{2}.results_,tags{1})
    rbsc = sc.resultbags{1};
else
    rbsc = sc.resultbags{2};
end

% set parameters
fieldBags = cell(size(tags));
md = java.security.MessageDigest.getInstance('MD5');
ws = weights; if iscell(ws), ws = cat(2,ws{:}); end
hash = md.digest([ws./max(ws(:)) double([tags{:}])]);
bi = java.math.BigInteger(1, hash);
md5 = ['s' char(bi.toString(35))]; md5 = md5(1:2:end);
if isempty(outName)
    wDir = [pDir a md5];
    resultDir = [wDir a 'project_results'];
    outName = [resultDir a 'fieldbag.jcm'];
    if exist(outName,'file'), out = {[];outName}; return; end
else
    b = strfind(outName,a);
    resultDir = outName;
    wDir = outName(1:b(end)-1);
end
if ~exist(resultDir, 'dir'), mkdir(resultDir); end

projectDir = sc.projectDir;
material = [projectDir a 'materials.jcmt'];
source = [projectDir a 'sources.jcmt'];
keys_t = sc.keys; keys_t.omega = omega;
if isfield(keys_t, 'relPerm') && isa(keys_t.relPerm, 'function_handle')
    keys_t.relPermittivity = keys_t.relPerm(keys_t.omega);
end
if contains(outName,{'background'})
    outfile = [wDir a 'sources.jcm'];
    jcmwave_jcmt2jcm(source,keys_t,'outputfile',outfile);
end
outfile = [wDir a 'materials.jcm'];
if exist(material,'file')
    jcmwave_jcmt2jcm(material,keys_t,'outputfile',outfile);
else
    copyfile(material(1:end-1),outfile);
end

for it = 1:length(tags)
    field = tags{it};
    if ~exist(field,'file')
        field = rbsc.results_.(tags{it}).result{1}.file;
    end
    ws = weights(it); if iscell(ws), ws = ws{1}; end
    ws = ['[' sprintf('(%.8e,%.8e),',[real(ws(:)) imag(ws(:))].')];
    ws(end) = ']';
    fieldBags{it} = sprintf(sc.fieldBag, field, ws);
end
fieldBags = [fieldBags{:}];

% create jcmpt file
jcmpt = [pDir a 'superposition.jcmpt'];
if ~exist(jcmpt,'file')
    superposition = sc.superposition;
    id = fopen(jcmpt, 'w');
    fprintf(id, superposition);
    fclose(id);
end
jcmpt = jcmpt(1:end-1);

% resultbag for superposition post process
if ~isfield(sc.quantities.Superposition,'resultbag')
    filename = [sc.resultbagDir a 'superposition.mat'];
    sc.quantities.Superposition.resultbag = jcmwave_resultbag(filename); 
end
resultbag = sc.quantities.Superposition.resultbag;

keys.fieldbags = fieldBags; keys.omega = omega;
keys.oPath = outName; tag = resultbag.get_tag(keys);
out = {jcmwave_solve(jcmpt,keys,resultbag,'workingdir',wDir);outName};
end
