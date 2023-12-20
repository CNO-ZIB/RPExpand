function updatePlots(RP,qs)
%UPDATEPLOTS update open plots
if nargin==1
    qs = {'ComplexPlane'}; 
elseif RP.up2date(5) || isempty(qs)
    return; 
else
    RP.up2date(5) = true;
end

plothandles = struct2cell(RP.plotSettings.fignumber);
for h = [plothandles{:}]
    if ~ishghandle(h), continue; end
    qiw = get(h,'UserData');
    if ~iscell(qiw) || isempty(qiw) || ~ischar(qiw{1})
        continue
    elseif ismember(qiw{1},qs)
        switch qiw{1}
            case RP.availableQuantities
                if ~all(inside(RP,qiw{end}))
                    qiw{end} = RP.expansionPoints; 
                end
                RP.plot(qiw{1:4})
            case 'Error'
                RP.plot(qiw{1:2});
            case 'ComplexPlane'
                RP.plot('ComplexPlane');
        end
    end
end
end

