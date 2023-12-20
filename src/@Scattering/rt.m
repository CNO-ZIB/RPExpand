function out = rt(quantity,tags)
%RT evaluates reflection and transmission coefficient, respectively
%   The evaluation is based on results from a Fourier transform. The
%   quantity is expected to have a field 
%   see also customizeQuantity

rb = quantity.resultbag; out = cell(size(tags)); 
header = rb.results_.(tags{1}).result{1}.header;
eps = header.RelPermittivity*quantity.parent.eps0;
mu = header.RelPermeability*quantity.parent.mu0;
nd = header.NormalDirection;
nd = ismember('XYZ',nd(end))*(1-2*(nd(1)=='m'));
c0 = 0.5*sqrt(eps/mu)/quantity.nm;
nOrders = length(rb.results_.(tags{1}).result{1}.N1);
contour = isfield(rb.results_.(tags{1}).keys,'circfield');
for it = 1:length(tags)
    k = rb.results_.(tags{it}).result{1}.K;
    res = rb.results_.(tags{it}).result{1};
    if nOrders~=length(res.N1) && contour
        error('The number of defraction orders must not change')
    end
    cf = res.ElectricFieldStrength{1};
    keys = rb.results_.(tags{it}).keys;
    if isfield(keys,'circfield')
        f1 = keys.circfield;
        keys.circfield = keys.field; keys.field = f1;
        if isfield(keys,'fieldIds')
            keys.fieldIds = keys.fieldIds(end:-1:1);
        end
        res_ = rb.get_result(keys); res_ = res_{1};
        cf_ = res_.ElectricFieldStrength{1};
    else, cf_ = cf; res_ = res;
    end
    c = c0*(k*nd(:))./sqrt(sum(k.*conj(res_.K),2));
    out{it} = c.'*sum(cf.*conj(cf_),2);
    if ~quantity.parent.derivatives, continue; end
    if isempty(quantity.parent.derivativeNames)
        fnames = fieldnames(res);
        fnames = fnames(startsWith(fnames,'d_'));
        quantity.parent.derivativeNames = ['FourierTransform';fnames];
    end
    fnames = quantity.parent.derivativeNames(2:end);
    ds = zeros(size(fnames));
    for it1 = 1:length(fnames)
        d = res.(fnames{it1}).ElectricFieldStrength{1};
        d_ = res_.(fnames{it1}).ElectricFieldStrength{1};
        ds(it1) = c.'*sum(d.*conj(cf_)+cf.*conj(d_),2);
    end
    out{it} = [out{it};ds.'];
end
out = cat(2,out{:});
end