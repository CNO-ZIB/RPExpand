function out = postprocess(sc,tags,qn)
%POSTPROCESS Apply post processes to the fields
%   out = POSTPROCESS(sc,tags,qn) A post process of JCMsuite is 
%   applied to each element of tags. The input parameter qn must evaluate 
%   to a valid fieldname of the property 'quantities'. The corresponding
%   quantity has a field 'jcmpt' containing the path to a post process.
%   The results from the post process are assigned tags again and a cell
%   array of tags, which belong to one contour, is passed to the function
%   handle provided with the field 'evaluate' of the quantity. This
%   function must take the quantity and the tags as arguments and return a
%   double array of the shape (d,n) where d is a custom data dimension and
%   n is the number of tags. The output is a cellaray that contains m
%   elements corresponding to m contours. In most cases m will be one.
%
%   If a Python expression is used for integration, the template file of
%   the post process must, within the section Python {...}, include the
%   definition IntegrationOrder = %(integrationOrder)e. This parameter
%   will be set automatically to twice the finite element degree of the
%   FEM solutions. In this section, you can also define the name of the
%   output quantity, e.g., IntegralName = "ElectricFieldEnergy". A Python
%   expression is currently always needed if you want to expand a quadratic 
%   functional like the energy of the electric field. As complex 
%   conjugation, being not holomorphic does not allow for a straightforward
%   application of this method, <a href="matlab: web(...
%   'https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.035432',...
%   '-browser')">circ fields</a> are used. Instead of
%   dot(E,E^*), dot(E,E^o) is computed, with E^o(w) = E(-w). To invoke the
%   field bag of a circ field within a Python expression you can define
%   FieldBagPath = "%(circfield)e". 
%   Please also refer to the <a href="matlab:
%   web(['https://docs.jcmwave.com/JCMsuite/html/ParameterReference' ...
%   '/d904047e6b00d14da719fb37a669e7eb.html'],'-browser')" 
%   >parameter reference</a> of JCMsuite and the 
%   documentation of the <a href="matlab:web(['https://docs.jcmwave.' ...
%   'com/JCMsuite/html/MatlabInterface/index.html'],'-browser')" 
%   >matlab interface</a> if you are interested in the 
%   definitions of the post processes in JCMsuite. The corresponding
%   project files must be stored in a directory whose path is passed to
%   the constructor of this class. Its default name is 'postprocesses' and
%   it is expected to be located in the same folder as @Scattering. Please
%   also refer to the existing post processes.
%
%   See also RieszProjection, Scattering,
%            jcmwave_resultbag, jcmwave_solve

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

% Check if a resource has been set up for computation
if isempty(jcmwave_daemon_resource_info)
    error('Please start a daemon');
end

% Make sure that the input parameters are as expected
if isnumeric(tags{1})
    out = sc.contributions(tags{1},qn); return
elseif ~isa(tags, 'cell')
    error('Cell expected but got %s', class(tags));
end
if isempty(qn)
    error('You must provide the quantity to be evaluated')
end

% For a dipole source the point evaluation for the emission already exists
if ismember(qn,{'DipoleEmission' 'PointEvaluation'})
    tags = tags(:,1); out = cell(1,length(tags)); q = sc.quantities.(qn);
    for it = 1:length(tags)
        out{it} = extractResult(q,tags{it},it==length(tags));
    end
    if iscell(out{end})
        sc.derivativeNames = out{end}{end}; out{end} = out{end}{1};
    end
    return;
end

if startsWith(qn,'reference')
    q = sc.quantities.(qn(10:end));
    if isfield(q,'reference')
        q = q.reference; 
    end
else
    q = sc.quantities.(qn);
end
keys = q.keys;
pP = q.jcmpt;
[~,pPName] = fileparts(pP);

% Get resultbags
if ~isfield(sc.pprbs,pPName)
    filename = [sc.resultbagDir filesep [pPName '.mat']];
    sc.pprbs.(pPName) = jcmwave_resultbag(filename);
end
rb = sc.pprbs.(pPName);
q.resultbag = rb;
    
jcmpt = q.jcmpt(1:end-1); % project file

nSets = size(tags,1); cField = size(tags,2)==2;
ids = cell(1,nSets); 
fieldIds = any(keys.fieldIds~=1);
if fieldIds, fIds = keys.fieldIds; end
for it1 = 1:nSets
    [files,w] = sc.tag2path(tags{it1,1});
    if cField, files_ = sc.tag2path(tags{it1,2}); end
    if isempty(w)
        wn = sc.eigenvalues(keys.fieldIds); w = wn.'; 
    end
    
    % Submit jobs to jcmwave_solve
    ids{it1} = zeros(size(files), 'int32');
    for it2 = 1:length(files)
        keys.integrationOrder = 2*sc.keys.finiteElementDegree+1;
        keys.oPath = 'out.jcm';
        keys.field = files{it2};
        if isfield(keys,'omega')
            keys.omega = w(it2);
        elseif isfield(keys,'lambda0')
            keys.lambda0 = 2*pi*sc.c0/w(it2);
        end
        if cField
            keys.circfield = files_{it2};
            if fieldIds, keys.fieldIds = [1 fIds(it2)]; end
            if fieldIds && all(tags{it1,1}{it2}(end-3:end)=='.jcm')
                keys.fieldIds = [fIds(it2) 1];
            end
        elseif fieldIds
            keys.fieldIds = 1;
        end
        ids{it1}(it2) = jcmwave_solve(jcmpt,keys,rb,'temporary','yes');
        tags{it1,1}{it2} = rb.get_tag(keys);
    end
end

% Collect results
jcmwave_daemon_wait([ids{:}],rb);
out = cell(size(ids)); 
args = cell(1,1+(nargin(q.evaluate)==3));
for it1 = 1:length(out)
    args{1} = tags{it1,1}; 
    if length(args)>1, args{end}=it1==length(out); end
    try out{it1} = q.evaluate(q,args{:});
    catch err
        fprintf(1,rb.results_.(tags{1}{1}).log.Log.Out);
        fprintf(1,'\n');
        rethrow(err);
    end
end

