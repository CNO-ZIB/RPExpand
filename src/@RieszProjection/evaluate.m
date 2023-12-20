function evaluate(RP)
%EVALUATE Evaluate scattering problems and update derived quantities
%   EVALUATE(RP) The scattering poblems along the contours are solved. The
%   fields of the property 'derivedQuantities' are updated.
%
%   see also RieszProjection/updateDerivedQuantities

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: April-2023

if isempty(RP.contours)
    if RP.useEigenvectors && ~RP.nPointsB
    else
        error('No contours: Please call the method ''getContours''.')
    end
end
d = max(0,length(RP.fields)-length(RP.contours)); fV = [RP.fields{:}];
if isempty(RP.undefinedValues)
    if ~isempty(RP.fields) || isempty(RP.contours)
        RP.updateDerivedQuantities; return; 
    end
    newV = RP.f(RP.contours); 
elseif d
    RP.fields(1:d) = [];
end

if any([RP.undefinedValues{:}])
    locs = find(cellfun(@any, RP.undefinedValues));
    frequencies = cell(1,length(locs));
    for it = 1:length(locs)
        l = locs(it);
        frequencies{it} = RP.contours{l}(RP.undefinedValues{l});
    end
    newV = RP.f(frequencies);
end

for it = 1:length(RP.contours)
    if isempty(RP.undefinedValues), RP.fields = newV; break; end
    locV = RP.definedValues{it};
    ndef = RP.undefinedValues{it};
    getEmpty = [class(RP.fields{1}) '.empty(0,length(ndef))'];
    RP.fields{it} = eval(getEmpty);
    if any(ndef)
        locnV = locs==it;
        RP.fields{it}(ndef) = newV{locnV};
    end
    if ~isempty(locV)
        RP.fields{it}(~ndef) = fV(locV);
    end
end

% update and reset quantities
RP.contours_f = RP.contours;
c = RP.contributions; ns = fieldnames(c).';
for n = ns, c.(n{1}) = rmfield(c.(n{1}),'residues'); end
RP.contributions = c;
RP.updateDerivedQuantities(d);
RP.undefinedValues = {};
RP.definedValues = {};
RP.up2date([1 5]) = [true false]; 
end

