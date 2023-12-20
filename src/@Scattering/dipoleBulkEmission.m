function varargout = dipoleBulkEmission(sc,omega,bE)
%DIPOLEBULKEMISSION evaluate the dipole bulk emission
%   emission = DIPOLEBULKEMISSION(sc,omega) The bulk emission at the given
%   frequencies omega is returned if the definition of the project
%   associated with the instance sc of the class Scattering contains a 
%   dipole emitter in a nonabsorbing material.
%   Known values are saved in the property 'bulkEmission' togeter with the
%   corrresponding frequencies

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

if ~isreal(omega)
    varargout{1} = getBulkEmission(sc,omega); return
end
bulkEmission = sc.quantities.DipoleEmission.reference.bulkEmission;

omega = sort(omega);
[idx,loc] = ismembertol(omega,bulkEmission(:,1));
if all(idx)
    if nargout, varargout{1} = bulkEmission(loc(idx),2); end
    return;
end

if nargin==3
    emission = bE;
else
    emission = zeros(size(omega));
    emission(idx) = bulkEmission(loc(idx),2);
    w0 = omega(~idx);
    emission(~idx) = getBulkEmission(sc,w0);
end
bulkEmission(loc(idx),:) = [];

bulkEmission = sortrows([[omega(:) emission(:)];bulkEmission]);
sc.quantities.DipoleEmission.reference.bulkEmission = bulkEmission;
if nargout, varargout{1} = emission; end
end

function res = getBulkEmission(sc,w0)
keys = sc.keys;
keys.relPermittivity = 1;
if (isfield(keys,'relPerm') && isa(keys.relPerm,'function_handle'))
    warning(['For the evaluation of the bulk emission with JCMsuite' ...
        ' the dispersion must be defined in the file materials.jcm']);
    keys = rmfield(keys,'relPerm');
end

projectDir = [sc.workingDir filesep 'bulk_emission'];
newFiles = false;
if ~exist(projectDir,'dir'), mkdir(projectDir); newFiles = true; end
if isfield(sc.quantities.DipoleEmission.reference,'rbblk')
    rb = sc.quantities.DipoleEmission.reference.rbblk;
else
    rb = jcmwave_resultbag([projectDir filesep 'resultbag.mat']);
    sc.quantities.DipoleEmission.reference.rbblk = rb; 
    newFiles = true; 
end

if newFiles
    materials = [sc.projectDir filesep 'materials.jcmt'];
    sources = [sc.projectDir filesep 'sources.jcmt'];
    layout = [sc.projectDir filesep 'layout.jcmt'];
    outfile = [projectDir filesep 'materials.jcm'];
    if exist(materials,'file')
        jcmwave_jcmt2jcm(materials,keys,'outputfile',outfile);
    else
        copyfile(materials(1:end-1),outfile);
    end
    outfile = [projectDir filesep 'layout.jcm'];
    if exist(layout,'file')
        jcmwave_jcmt2jcm(layout,keys,'outputfile',outfile);
    else
        copyfile(layout(1:end-1),outfile);
    end
    outfile = [projectDir filesep 'sources.jcm'];
    keys.omega = w0(1);
    jcmwave_jcmt2jcm(sources,keys,'outputfile',outfile);
    jcmwave_geo(projectDir);
end
projectFile = [projectDir filesep 'project.jcmpt'];
if ~exist(projectFile,'file')
    id = fopen(projectFile,'w');
    fprintf(id,sc.pP,projectDir);
    fclose(id);
end

projectFile = projectFile(1:end-1);
keys.omega = w0;
try
    id = jcmwave_solve(projectFile,keys,rb,'temporary','yes');
catch me
    if ~contains(me.message,'resultbag is invalidated'), rethrow(me); end
    rb.reset; id = jcmwave_solve(projectFile,keys,rb,'temporary','yes');
end
jcmwave_daemon_wait(id,rb);
res = rb.get_result(keys);
try res = res{1}.BulkEmission{1};
catch
    logs = rb.get_log(keys);
    error(logs.Log.Error);
end
end

