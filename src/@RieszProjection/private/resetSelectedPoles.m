function resetSelectedPoles(RP,selectedPoles)
%RESETSELECTEDPOLES set selected poles, autoselect if required
%   If the property autoselect is true, all poles inside the contour are
%   selected automatically

if isempty(selectedPoles)
    selectedPoles = RP.selectedPoles;
end
inC = {}; n_old = length(RP.selectedPoles);
if RP.autoselect && ~isempty(RP.contours)
    inC = find(inside(RP,RP.poles));
    sel = true(1,length(selectedPoles));
    for it = 1:length(selectedPoles)
        sel(it) = ~all(ismember(selectedPoles{it},inC));
    end
    selectedPoles(sel) = [];
    inC = num2cell(inC(~ismember(inC,[selectedPoles{:}])));
end
selectedPoles = [selectedPoles inC];
if isequal(selectedPoles, RP.selectedPoles), return; end
n = max(length(selectedPoles),1);
if n~=length(RP.nPoints_)
    if any(RP.nPoints_~=RP.nPoints_(1))
        warning(['The value of property '...
            '''nPoints'' has been changed.'])
    end
    RP.nPoints_ = repmat(min(RP.nPoints),1,n);
end
if n~=length(RP.radius_)
    if any(RP.radius_~=RP.radius_(1))
        warning(['The value of property '...
            '''radius'' has been changed.'])
    end
    RP.radius_ = repmat(min(RP.radius),1,n);
end
RP.selectedPoles_ = selectedPoles;
RP.nModes = length(selectedPoles);
RP.contributions = struct;
if RP.nPoints(1) && n_old~=n, RP.getContours(true);
elseif isfield(RP.plotSettings,'fignumber') && n_old~=n
    if isfield(RP.plotSettings,'index')
        RP.plotSettings = rmfield(RP.plotSettings,'index');
    end
    if isfield(RP.plotSettings.fignumber,'ComplexPlane') && ...
            ishandle(RP.plotSettings.fignumber.ComplexPlane)
        RP.plot('ComplexPlane')
    end
end
end