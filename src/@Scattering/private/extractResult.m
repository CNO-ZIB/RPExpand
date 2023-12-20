function out = extractResult(quantity,tags,last)
%EXTRACTRESULT extract results from post processes
%   If you want to add a custom quantity to the interface, you must specify
%   how the results of the post process is used to evaluate the target
%   quantity. For this purpose each quantity that is added to the interface
%   must be equipped with a function handle that takes the name of the 
%   quantity and tags, which refer to the resultbag, and returns
%   the corresponding values of the target quantity. You can either define
%   this function sperately or add a suitable case to the switch statement
%   below. This function is used as the default if no function handle is
%   passed to the method <a href="matlab: 
%   help('Scattering/quantity')">quantity</a> explicitly. 

switch quantity.id
    case {'referenceFarFieldIntegral' 'FarFieldIntegral'}
        out = integrateFarfield(quantity,tags,last);
    case {'referenceRadiationPattern' 'RadiationPattern'}
        out = pointingVector(quantity,tags);
    case 'referenceDipoleEmission'
        out = referenceDipoleEmission(quantity,tags);
    case {'DipoleEmission','PointEvaluation'}
        out = dipoleEmission(quantity,tags,last);
    case 'ModeCoupling'
        out = coupling(quantity,tags);
    case 'FourierTransform'
        out = fourierTransform(quantity,tags);
    case 'CartesianExport'
        out = cartesianExport(quantity,tags);
    otherwise
        out = getResult(quantity,tags);
end
end

function out = getResult(quantity,tags)
name = quantity.id;
if startsWith(name,'reference'), name = name(10:end); end
out = cell(size(tags));
for it = 1:length(tags)
    v = quantity.resultbag.results_.(tags{it}).result{end};
    out{it} = v.(name){end}(:);
    if any(contains(fieldnames(v),'DomainId'))
        if isfield(v,'DomainIds') && isfield(quantity,'DomainIds')
            out{it} = out{it}(ismember(v.DomainIds,quantity.DmainIds));
        end
        out{it} = sum(out{it});
    end
end
out = cat(2,out{:});
end

function out = dipoleEmission(quantity,tags,last)
resultbag = quantity.resultbag;
keys = resultbag.results_.(tags{1}).keys;
out = cell(size(tags)); 
for it = 1:length(tags)
    result = resultbag.results_.(tags{it}).result;
    fname = [result{1}.file(1:end-12) 'points.jcm'];
    if last && it==length(tags) && quantity.parent.derivatives
        [fvs,quantity.parent.derivativeNames] = Scattering.load(fname);
    elseif quantity.parent.derivatives
        fvs = Scattering.load(fname);
    else
        fvs = result{end}.field{1}.';
    end
    if isfield(keys,'strength') && strcmp(quantity.id,'DipoleEmission')
        out{it} = -0.5*sum(fvs.*conj(keys.strength(:)).',2);
    else
        out{it} = sum(fvs.*quantity.component(:).',2);
    end
end
out = cat(2,out{:}); 
end

function out = referenceDipoleEmission(quantity,tags)
out = zeros(size(tags)); bE = zeros(length(tags),2);
for it = 1:length(tags)
    v1 = quantity.resultbag.results_.(tags{it}).result{end};
    out(it) = v1.DipoleEmission;
    bE(it,2) = v1.DipoleBulkEmission;
    bE(it,1) = quantity.resultbag.results_.(tags{it}).keys.omega;
end
dipoleBulkEmission(quantity.parent,bE(:,1),bE(:,2));
end

function out = integrateFarfield(quantity,tags,last)
resultbag = quantity.resultbag; out = cell(size(tags));
persistent results resize;
if isempty(resize), resize = false; results = resultbag.results_; end
for it = 1:length(tags)
    tag = tags{it};
    farfield = resultbag.results_.(tag).result;
    % Check if the integration has already been done
    if isnumeric(farfield)
        out{it} = farfield;
    else
        keys_ = resultbag.results_.(tag).keys;
        pwt = RieszProjection.weights15;
        WT = [pwt(end:-1:2), pwt];
        r = keys_.radius; ndx = keys_.ndx-1;
        theta = reshape(keys_.gridPointsTheta,15,[]);
        sz = [size(theta) length(keys_.gridPointsPhi) 3];
        field1 = reshape(farfield{1}.ElectricFieldStrength{1},sz);
        if isfield(keys_,'circfield')
            f1 = keys_.circfield;
            keys_.circfield = keys_.field; keys_.field = f1;
            if isfield(keys_,'fieldIds')
                keys_.fieldIds = keys_.fieldIds(end:-1:1);
            end
            field2 = resultbag.get_result(keys_);
            field2 = reshape(field2{1}.ElectricFieldStrength{1},sz);
        else, field2 = field1;
        end
        hlfhp = deg2rad(keys_.halfhp(1));
        hlfht = deg2rad(keys_.halfht); c1 = r^2*hlfhp;
        eps = farfield{1}.header.RelPermittivity*Scattering.eps0;
        mu = farfield{1}.header.RelPermeability*Scattering.mu0;
        S = 0.5*sqrt(eps/mu)*sum(field1.*conj(field2),4).*sind(theta);
        P = c1.*sum(S.*repmat(reshape(WT,1,1,[]),1,1,size(S,3)/15),3);
        out{it} = sum(P(:,1:ndx).*WT.'.*hlfht(1:ndx),[1 2]);
        if size(P,2)>ndx, out{it}(2,1) = sum(P.*WT.'.*hlfht,[1 2]); end
        resize = true;
    end
    results.(tag).result = out{it}; % replace fields whith integrals
end
out = cat(2,out{:});
if last % reload resultbag with smaller size
    if resize
        fprintf('\nPrepare resized resultbag.')
        results.fieldnames = resultbag.fieldnames_;
        results.keys = struct;
        results.source_files = resultbag.source_files_;
        delete(resultbag.filepath_)
        save(resultbag.filepath_,'-struct','results')
        quantity.resultbag = jcmwave_resultbag( ...
            resultbag.filepath_, resultbag.fieldnames_); fprintf(1,'\n');
        [~,pPName] = fileparts(quantity.jcmpt);
        quantity.parent.pprbs.(pPName) = quantity.resultbag;
    end
    results = []; resize = [];
end
end

function out = pointingVector(quantity,tags)
resultbag = quantity.resultbag; out = cell(size(tags));
for it = 1:length(tags)
    tag = tags{it};
    keys_ = resultbag.results_.(tag).keys;
    farfield = resultbag.results_.(tag).result;
    field1 = farfield{1}.ElectricFieldStrength{1};
    if isfield(keys_,'circfield')
        f1 = keys_.circfield;
        keys_.circfield = keys_.field; keys_.field = f1;
        if isfield(keys_,'fieldIds')
            keys_.fieldIds = keys_.fieldIds(end:-1:1);
        end
        field2 = resultbag.get_result(keys_);
        field2 = field2{1}.ElectricFieldStrength{1};
    else, field2 = field1;
    end
    eps = farfield{1}.header.RelPermittivity*Scattering.eps0;
    mu = farfield{1}.header.RelPermeability*Scattering.mu0;
    out{it} = 0.5*sqrt(eps/mu)*sum(field1.*conj(field2),2);
    out{it} = out{it}(:);
end
out = cat(2,out{:});
end

function out = fourierTransform(quantity,tags)
resultbag = quantity.resultbag; out = cell(size(tags)); 
ndx = 1; cmp = [1 1 1]; % defaults
if isfield(quantity,'diffractionOrder'), ndx = 0; end
if isfield(quantity,'component'), cmp = quantity.component; end
dOrder = quantity.diffractionOrder;
if isscalar(dOrder), dOrder = dOrder([1 1]); end
nOrders = length(resultbag.results_.(tags{1}).result{1}.N1);
contour = isfield(resultbag.results_.(tags{1}).keys,'circfield');
for it = 1:length(tags)
    tag = tags{it};
    res = resultbag.results_.(tag).result{1};
    if nOrders~=length(res.N1) && contour
        error('The number of defraction orders must not change')
    end
    if ~ndx, ndx = find(res.N1==dOrder(1)&res.N2==dOrder(2),1); end
    out{it} = res.ElectricFieldStrength{1}(ndx,:);
    out{it} = out{it}*cmp(:);
    if ~quantity.parent.derivatives, continue; end
    if isempty(quantity.parent.derivativeNames)
        fnames = fieldnames(res);
        fnames = fnames(startsWith(fnames,'d_'));
        quantity.parent.derivativeNames = ['FourierTransform';fnames];
    end
    fnames = quantity.parent.derivativeNames(2:end);
    ds = zeros(size(fnames));
    for it1 = 1:length(fnames)
        d = res.(fnames{it1}).ElectricFieldStrength{1}(ndx(1),:);
        ds(it1) = d*cmp(:);
    end
    out{it} = [out{it};ds.'];
end
out = cat(2,out{:});
end

function out = coupling(quantity,tags)
resultbag = quantity.resultbag; out = cell(size(tags));
for it = 1:length(tags)
    tag = tags{it}; keys_ = resultbag.results_.(tag).keys;
    ndx = keys_.fieldIds; result = resultbag.results_.(tag).result;
    c1 = result{3}.ModeCoefficients{ndx(1)};
    if isfield(keys_,'circfield')
        if isfield(keys_,'lambda0')
            keys_.lambda0 = conj(keys_.lambda0); 
        end
        f1 = keys_.circfield;
        keys_.circfield = keys_.field; keys_.field = f1;
        if isfield(keys_,'fieldIds')
            keys_.fieldIds = keys_.fieldIds(end:-1:1);
        end
        result_ = resultbag.get_result(keys_);
        c2 = result_{3}.ModeCoefficients{ndx(end)};
        E1 = real([result{2}.field{:}]); E2 = real([result_{2}.field{:}]);
        u1 = complex(E1(1,1:2),E1(2,1:2)).'; 
        u2 = complex(E2(1,1:2),E2(2,1:2)).';
        phi = diff(angle([u2 u1]),1,2);
        if length(c2)==2
            tm = eye(2);
            if abs(abs(diff(phi))-pi)<1e-6
                [~,ndx] = min(abs(phi)); phi = phi(ndx);
                d = ones(2,1)*-1; d(ndx) = 1; tm = diag(d);
            else
                phi = phi(1);
            end
            tm = tm*[cos(phi) -sin(phi); sin(phi) cos(phi)];
            c2 = tm*c2;
        elseif all(abs(abs(phi)-pi)<1e-3)
            c2 = -c2;
        end
    else, c2 = c1;
    end
    out{it} = sum(c1.*conj(c2));
end
out = cat(2,out{:});
end

function out = cartesianExport(quantity,tags)
resultbag = quantity.resultbag; out = cell(size(tags));
for it = 1:length(tags)
    result = resultbag.results_.(tags{it}).result{1};
    field = result.field{1};
    out{it} = field(:);
end
out = cat(2,out{:});
end