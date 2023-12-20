function out = computeContributions(RP, varargin)
%COMPUTECONTRIBUTIONS compute the individual contributions
%   You will usually not need to use this method directly but rather use
%   computeExpansion, plot, computePoles or computeZeros
%
%   see also computeExpansion, plot, computePoles, computeZeros

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: April-2023

persistent parser;
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/integrate';
    addRequired(parser,'quantity',@(x)validateattributes(x,{'char'},{}));    
    addOptional(parser,'w0',0,...
        @(x)validateattributes(x,{'numeric'},{}));
    addParameter(parser,'precision',0,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative'}));  
    addParameter(parser,'upperbound',128,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'})); 
    addParameter(parser,'index',0,...
        @(x)validateattributes(x+1,{'numeric'},{'nonnegative','scalar'})); 
end
parse(parser,varargin{:});
vars = parser.Results; q = vars.quantity;
if ~vars.precision, vars.precision = RP.precision; end % for adaptivity

% if necessary update eigenvalues
if ~(RP.up2date(2)||RP.up2date(3)), RP.computePoles(q); end
% if necessary update scattering solutions and derived quantities
if ~RP.up2date(1), RP.evaluate; end
if isfield(RP.contributions,q) && isfield(RP.contributions.(q),'residues')
    vars.index = -1; 
elseif RP.useEigenvectors
    RP.contributions.(q) = struct; RP.updateDerivedQuantities;
end

% if only the background is required check if that should be done using
% interpolation 
if isfield(RP.contributions,q) && ~RP.nPointsB 
    if ~isfield(RP.contributions.(q),'interpolation')
        quantity = RP.quantities.(q); q2 = quantity.quadratic;
        sz = [max(2,RP.nInterpolationNodes) quantity.m];
        residues = RP.contributions.(q).residues;
        residues = reshape(residues,1,size(residues,2),[],1+q2);
        wn = reshape(RP.poles([RP.selectedPoles{:}]),1,1,[]);
        if q2, wn = cat(4,wn,conj(wn)); end
        w0 = RP.expansionPoints; 
        w0 = linspace(w0(1),w0(end),max(sz(1),2)).';
        if RP.nInterpolationNodes
            ref = RP.computeReference(q,w0);
            ms = cat(3,sum(-1/(2i*pi)*residues./(wn-w0),4),zeros(sz));
        else
            ms = zeros([sz size(residues,3)+1]);
        end
        if isfield(quantity,'evaluate')
            out = {ms w0};
            ms = quantity.evaluate(out{1:nargin(quantity.evaluate)});
        end
        if ~RP.nInterpolationNodes, ref = zeros(size(ms,1,2)); end
        coefficients = spline(w0,ref-real(sum(ms,3)));
        RP.contributions.(q).interpolation = coefficients;
    end
    out = ppval(RP.contributions.(q).interpolation,vars.w0);
    return;
end
vars.quadgk = RP.shape{3}>1; % number of gk nodes

nContours = length(RP.contours);
q2 = RP.quantities.(q).quadratic;
if ~isempty(RP.poles)
    poles = [{} cellfun(@(x){RP.poles(x)},RP.selectedPoles)];
    if q2, poles = [poles cellfun(@(x){conj(x)},poles)]; end
    nPoles = cumsum(cellfun(@length,poles));
    if isempty(nPoles), nPoles = 0; end
    if nContours==1, poles = {[poles{:}]}; nPoles = nPoles(end); end
else
    nPoles = 0;
end
nContributions = nPoles(end)+1;
c = 1; p = 1; 
if vars.index % the contours to be integrated
    c = vars.index; inds = vars.index;
    if c<0, c = nContributions; inds = nContributions;
    elseif c<length(nPoles) && ~strcmpi(q, 'field')
        inds = nPoles(c):nPoles(c+1)-1;
    end
    if nContours==1, c=1; p = vars.index; end
else
    inds = 1:nContributions;
end

res = cell(size(inds));
for it1 = inds
    vars.p = p;
    if it1==nContributions
        c = nContours; ps = [];
        if ~RP.nPointsB, continue; end
    else
        ps = poles{c}; p = p+1;
    end
    res{min(it1,end)} = RP.quad(ps,c,vars); % integrate
    if it1==nPoles(min(c,end)), c = min(c+1,nContours); p = 1; end
end
out = 1/(2i*pi)*res{end}; % The background contribution
if ~vars.index && ~isempty(RP.selectedPoles)
    RP.contributions.(q).residues = cat(3,res{1:end-1});
elseif isempty(RP.selectedPoles)
    RP.contributions.(q).residues = zeros(1,RP.quantities.(q).m,0);
end
if any(isinf(vars.precision)), return; end

% adaptive refinement
% check if the target precision is reached 
tol = vars.precision; % tol(1) absolute and tol(end) relative tolerance
if isscalar(tol), tol = [tol tol]; end 
indices = [1:min(RP.nModes,nContours) nContours+(0:1)]; % contours and ref
finished = RP.computeError(q,indices,tol);
if isfield(RP.reference,q)
    if finished(end), return; end % comparison with physical quanity
    finished = finished(1:end-1); gl = false;
end
if isempty(finished) || (all(finished) && length(finished)==nContours-1)
    % add Kronrod nodes to background, which was integrated with Gauss
    finished = [finished false]; gl = true; % global refinement
elseif length(finished)==nContours && RP.nPoints_(1)
    finished = [finished true]; % refine modal contours first
elseif all(finished), return;
end
% edit number of points if required
nPnts = RP.nPoints_; 
if nPnts(1)
    nPnts(~finished(1:RP.nModes)) = nPnts(~finished(1:RP.nModes))*2;
    finished(1:RP.nModes) = finished(1:RP.nModes)|nPnts>vars.upperbound;
    RP.nPoints(~finished(1:RP.nModes)) = nPnts(~finished(1:RP.nModes));
    if all(finished), return; end
end

if length(finished)>RP.nModes || ~RP.nPoints_(1)
    if vars.quadgk && any(~finished(RP.nModes+1:end))
        if ~gl
            v = RP.shape{1}; cr = []; sp = RP.shape{2};
            if ismember(sp,'ec'), cr = v(1:2); v = v(3:end); end
            ndx = ~finished(RP.nModes+1:end); 
            if ~any(ndx), return; end
            k = 1; c = RP.includeConjugatedPoles;
            if ismember(sp,'ec') && c, ndx = ndx|ndx(end:-1:1); end
            c = c && sp=='p';
            v_new = zeros(1,length(v)+sum(ndx)+c*sum(ndx)-1);
            for it = 1:length(ndx)
                v_new(k) = v(it); k = k+1;
                if ndx(it)
                    v_new(k) = (v(it+1)+v(it))/2; k = k+1;
                    if c, v_new(k) = conj(v_new(k-1)); k = k+1; end
                end
            end
            if c % sort the vertices to form a polygon
                phi = mod(2*pi+eps(10)+angle(v_new-mean(v_new)),2*pi);
                [~,idx] = sort(phi);
                d = abs(diff(v_new(idx)));
                idx = idx([d>eps(max(d))*10 true]);
                v_new = v_new(idx);
            end
            if c, endPoint = v_new(1); else, endPoint = v(end); end
            RP.shape{1} = [cr v_new endPoint];
        end
        RP.getContours('gk');
        RP.nPointsB_ = numel(RP.contours{end});
        RP.up2date(2) = RP.up2date(3) && ~RP.nPoints;
    elseif ~finished(end)
        RP.nPointsB=RP.nPointsB*2;
    end
end
if all([finished numel(RP.contours{1})>vars.upperbound])
    return;
end

out = RP.computeContributions(q, vars.w0, ...
    'precision', vars.precision, 'upperbound', vars.upperbound);
end
