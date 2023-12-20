function updateDerivedQuantities(RP, d)
%UPDATEDERIVEDQUANTITIES Update all quantities
%    UPDATEDERIVEDQUANTITIES(RP) All quantities saved as fields of the
%        property 'derivedQuantities' are updated.
%
%    UPDATEDERIVEDQUANTITIES(RP,d) If the number of contours has changed
%        the cell array is shrunk such that its length matches the number 
%        of modes.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

% The derived quantities may either be scalars or vectors of length m. They
% are evaluated at p*q integration points, where, in the case of a higer-order
% quadrature, q>1 is the number of subintervals. The dimension of the 
% derived quantity of a given contour must be (1,m,p,q). The first
% dimension is reserved for the number of expansion points. If derivatives
% with respect to design parameters are available, they must be stored in
% the data dimension.

if nargin==1, d=0; end % d elements will be removed to shrink the size
derivedQs = RP.derivedQuantities;
useEigenvectors = RP.useEigenvectors && d>-1;
if useEigenvectors, derivedQs = RP.contributions; end
if isempty(derivedQs), return; end
quantityNames = fieldnames(derivedQs); 
gk = ~isempty(RP.contours) && RP.shape{3}>1;
nContours = length(RP.contours_f); 
nMC = min(nContours,RP.nModes+1);

for it1 = 1:length(quantityNames)
    q = quantityNames{it1};
    getAll = isempty(derivedQs.(q));
    
    % get additional parameters
    isLinear = ~RP.quantities.(q).quadratic;
    
    % if normalized eigenvectors are used, get function handles for modal
    % contributions if not yet existing
    if useEigenvectors && ~isfield(derivedQs.(q),'residues')
        w_n = RP.poles([RP.selectedPoles{:}]);
        RP.contributions.(q).residues = RP.f({w_n},q);
        if ~isfield(RP.derivedQuantities,q)
            continue
        else
            getAll = isempty(RP.derivedQuantities.(q));
        end
    end
    
    % get existing results and iterate over the contours
    if ~getAll
        if isempty(RP.undefinedValues), continue; end
        undefined = RP.undefinedValues([1:nMC-1 end]);
        quantity = RP.derivedQuantities.(q);
        quantity = cellfun(@(x){reshape(x,size(x,2),[])},quantity);
        if isLinear, shift = sum(cellfun(@length,quantity(1:end-1))); end
        quantity = cat(2,quantity{:});
        if d>0, RP.derivedQuantities.(q)(1:round(d/(1+isLinear))) = []; end
    end
    fs = cell(nContours,1+RP.includeConjugatedPoles);
    for it2 = 1:nMC
        if it2==nMC, it = nContours; else, it = it2; end
        vsp = RP.fields{it};
        if isempty(vsp) || (useEigenvectors && it2<nMC), continue; end
        if getAll
            ndef = [];
        elseif any(undefined{it2})
            ndef = undefined{it2};
            vsp = vsp(ndef);
        else
            continue
        end
        if RP.includeConjugatedPoles
            if it2<nMC, vsn = reverseV(vsp,it,ndef); end
            vsp_ = RP.fields{mod(it2-1,nMC)+nMC};
            if ~isempty(ndef), vsp_ = vsp_(ndef); end
            vs = {vsp reverseV(vsp_,it,ndef)};
        else, vs = {vsp};
        end
        fs(it,:) = vs;
        
        % quadratic case
        if ~isLinear && it2<nMC, fs(it2+nMC-1,:) = {vsp_,vsn}; end
    end
    
    notEmpty = ~cellfun('isempty',fs(:,1)); v = cell(size(fs,1),1);
    v(notEmpty) = RP.f(fs(notEmpty,:),q);
    
    % distribute the results
    for it2 = 1:nContours
        if it2==nMC&&isLinear, it = nContours; else, it = it2; end
        if isLinear && it2>nMC, break; end
        sz = size(RP.contours_{it});
        if ~getAll
            ndef = RP.undefinedValues{it};
            locV = RP.definedValues{it};
            if (it2==nMC&&isLinear) && ~isempty(locV) ...
                    && RP.includeConjugatedPoles
                locV = locV-shift;
            end
            RP.derivedQuantities.(q){it2} = ...
                zeros([1,size(quantity,1),sz]);
            RP.derivedQuantities.(q){it2}(1,:,~ndef) = ...
                quantity(:,locV);
            if isempty(v{it2}), continue; end
        end
        if getAll 
            % axis 1: w0, axis 2: data, axis 3: nodes, axis 4: subs
            sz_ = [1 size(v{it2},1) sz];
            RP.derivedQuantities.(q){it2} = reshape(v{it2},sz_);
        else
            RP.derivedQuantities.(q){it2}(1,:,ndef) = v{it2};
        end
    end
end

    function vr = reverseV(values,ndx,ndef)
        c1 = RP.contours{ndx}(:);
        M = max(abs(c1));
        if ~isempty(ndef), c1 = c1(ndef); end
        c1 = c1([1 end]);
        vr = values(end:-1:1);
        if gk && ndx==length(RP.contours)
            ngk = RP.shape{3};
            if imag(c1(2))>0 && isempty(ndef)
                vr = circshift(vr,-ngk);
            elseif imag(c1(2))>0
                n = sum(mod(find(~ndef)+ngk/2,length(ndef))<ngk);
                vr = circshift(vr,-ngk+n);
            end
        else
            if ndx<nMC
                c2 = RP.contours{ndx+nMC-1}; 
                if ~isempty(ndef), c2 = c2(ndef); end
                c2 = c2(1);
            else
                c2 = c1(1);
            end
            if imag(c2)+imag(c1(1))<2*eps(M)
                vr = circshift(vr,1);
            end
        end
    end
end

