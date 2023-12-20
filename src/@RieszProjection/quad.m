function res = quad(RP,g_p,c,vars)
%QUAD apply quadrature rule. Either the trapezoidal rule or Gauss-Kronrod.
%   The input p_g must be the analytic contribution g to the integral or
%   poles. In the first case a pointwise multiplication with the weights
%   will be performed, i.e., weights.*g, in the second case an analytic
%   contributions will be evaluated based on the poles and the index to one
%   of the poles provided with vars.p.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

n = abs(c); err = 0; if isfield(vars,'err'), err = vars.err; end
if isfield(vars,'p'), p = vars.p; else, p = 0; end % index for poles
integrate = ~strcmp(vars.quantity,'weights');
if integrate, m = RP.quantities.(vars.quantity).m; end % data dimension
if isfield(vars,'dvs') && vars.dvs, m = vars.dvs+1; end % derivatives

w = reshape(RP.contours{n},[1 1 size(RP.contours{n})]);
if ~p
    g = g_p;
elseif isempty(g_p)
    g = 1./(w-vars.w0);
else
    idx_wn = zeros(size(g_p),'logical');
    idx_wn(p) = true;
    wn = g_p(idx_wn);
    wns = g_p(~idx_wn);
    if ~isempty(wns)
        g = prod((w-wns)./(wn-wns),2);
    else
        g = 1;
    end
end

% Dimensions of 'values':
% axis 1: w0, axis 2: data, axis 3: nodes, axis 4: subintervals
if strcmpi(vars.quantity, 'field')
    values = RP.fields{min(n,end)};
elseif isfield(RP.derivedQuantities,vars.quantity)
    values = RP.derivedQuantities.(vars.quantity){min(n,end)};
    values = values(:,1:m,:,:); % ignore derivatives if not required
elseif ~integrate
    values = NaN;
else
    error('No data available for integrating the %s',vars.quantity)
end
if c<0 % inverse quantity (only background contour)
    values(1,1,:) = 1./values(1,1,:);
    if isfield(vars,'dvs') && vars.dvs % derivative of quotient
        values(1,2:end,:) = -values(1,1,:).^2.*values(1,2:end,:);
    end
end

if isempty(values), res = NaN; return; end
weights = RP.weights{min(n,end)};

if iscell(weights)
    weights = reshape(weights{1+logical(err)},[1 1 size(weights{1})]).*g;
else
    weights = reshape(weights,[1 1 size(weights)]).*g;
    if err % return error instead of integral
        weights = repmat(weights,1,1,1,err);
        for it=1:err
            n = 2^it; % halve the number of integration points in each step
            weights(:,:,1:n:end,it) = -weights(:,:,1:n:end,it)*(n-1);
        end
    end
end

if ~integrate, res = weights; return; end

% integrate
if RP.shape{3}>1 && ~err
    res = sum(weights.*values,[3 4]);
else
    res = sum(weights.*values,3);
end
end
