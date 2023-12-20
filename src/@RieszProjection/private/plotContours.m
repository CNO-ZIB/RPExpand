% plot contours, expansion points and eigenvalues

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022
function fig = plotContours(RP, vars)
if vars.maxIndex > 1
    nModes=length(RP.selectedPoles);
else
    nModes = 0; % if there is not any contour yet or no poles are selected
end
% plot
if ~vars.fignumber, RP.plotSettings.fignumber.ComplexPlane = 1; end
fig = figure(RP.plotSettings.fignumber.ComplexPlane);
clf(vars.fignumber);
ax1 = axes(fig); ax1.NextPlot = 'add'; ax1.Box = 'on'; axis(ax1,'equal')
ax1.TickLabelInterpreter = vars.opts{2}; ax1.FontSize = vars.opts{4};
if strcmpi(vars.interpreter,'latex')
    xlabel(ax1,'Re($\omega$) $\left[\mathrm{s}^{-1}\right]$',vars.opts{:})
    ylabel(ax1,'Im($\omega$) [$\mathrm{s}^{-1}$]',vars.opts{:})
else
    xlabel(ax1,'Re(\omega) [s^{-1}]',vars.opts{:})
    ylabel(ax1,'Im(\omega) [s^{-1}]',vars.opts{:})
end
if vars.index
    inds = vars.index;
    if all(inds==nModes+1), inds = 1:nModes+1; end
else
    inds = 1:nModes+1;
end
if RP.includeConjugatedPoles
    inds = [inds(inds<nModes+1) nModes+inds(end)];
end
if ~isempty(RP.poles)
    ndx = ismember(1:length(RP.poles),[RP.selectedPoles{:}]); 
    args = {'MarkerSize',6,'LineWidth',1,'DisplayName','$\omega_n$'};
    plot(ax1,RP.poles(ndx), 'rx',args{:});
end
if ~isempty(RP.expansionPoints)
    ln = 'k-'; if length(RP.expansionPoints)<3, ln = 'k.-'; end
    l2 = plot(ax1,RP.expansionPoints, zeros(size(RP.expansionPoints)),...
        ln,'DisplayName','$\omega_0$','LineWidth',1);
    l2.MarkerSize = 10;
end
for it = inds
    contour = 0;
    try contour = RP.contours{it}; catch, break; end
    if isempty(contour), continue; end
    if it==inds(end)
        if isscalar(contour) && isreal(contour)
            plot(ax1,contour,0,'ko','DisplayName','Integration Points');
        else
            plot(ax1,contour(:),'k.','DisplayName','Integration Points');
            [v,s] = RP.shape{1:2};
            if s == 'c' || s == 'e' % circular shape
                c = v(1); r = v(2);
                alpha = linspace(0,2*pi,10*numel(contour));
                if s == 'c'
                    contour = c + r*exp(1i*alpha);
                elseif s == 'e'
                    E = RP.shape{4}; c = E(1)+1i*E(2);
                    if RP.includeConjugatedPoles, alpha = alpha-E(5); end
                    A = E(3); B = E(4); a = alpha+E(5); b = alpha-E(5);
                    contour = c + A*exp(1i*a)+B*exp(-1i*b);
                end
            else
                contour = [contour(:);contour(1)];
            end
            plot(ax1,contour,'k')
        end
    else
        plot(ax1,contour(:),'r.','DisplayName','');
    end
end
selection = ax1.Children(2-isreal(contour):end);
if length(selection)>3
    selection = selection([1 end-1:end]);
end
if ~isempty(RP.poles) && any(~ndx) % should not appear in the legend
    plot(ax1,RP.poles(~ndx),'x',args{1:4},'Color',[0.5 0.5 0.5]);
end
cindex = max(vars.index); 
if cindex>length(RP.selectedPoles), cindex = length(RP.contours); end
if any(vars.ylim)
    ylim(vars.ylim)
elseif vars.index
    M = max(imag(RP.contours{cindex})); 
    m = min(imag(RP.contours{cindex}));
    shift = (M-m)/10; ylim([m-shift M+shift]);
end
if any(vars.xlim)
    xlim(vars.xlim)
elseif vars.index
    shift = (M-m)/2;
    M = max(real(RP.contours{cindex}));
    m = min(real(RP.contours{cindex}));
    shift = max(shift-(M-m)/2,(M-m)/10); xlim([m-shift M+shift]);
end
legend(ax1,selection,'location',vars.legendLocation,vars.opts{:})
legend(ax1,'boxoff')
if isempty(RP.contours_), return; end

pos = get(ax1, 'Position');
if vars.contourNumbers && ~isempty(RP.selectedPoles)
    ax2 = axes('Position',pos,'XAxisLocation','top',...
        'Color','none','XColor','k','YColor','none');
    ax2.TickLabelInterpreter = vars.opts{2}; ax2.FontSize = vars.opts{4};
    axis(ax2,'equal')
    xlim(get(ax1,'xlim')); ylim(get(ax1,'ylim'));
    inds = 1:length(RP.selectedPoles);
    if all(RP.selectedPoles{1})
        xticklabels(arrayfun(@num2str,inds,'UniformOutput',false))
        ps = cellfun(@(x) x(1), RP.selectedPoles);
        xticks(ax2,real(RP.poles(sort(ps(inds)))));
        yticks(ax2,[]);
    end
    linkaxes([ax1 ax2],'xy'); 
    fig.CurrentAxes = ax1; ax1.Box = 'off';
end
end
