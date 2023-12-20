% plot the error of the contour integrations

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022
function fig = plotError(RP,vars)

if ~any(vars.index), vars.index = length(RP.contours)+(0:1); end
errortype = lower(vars.errortype(1));
if errortype == 'r', label_y = 'Maximal Relative Error';
elseif errortype == 'a', label_y = 'Maximal Absolute Error';
else, label_y = 'Error'; 
end
errortype = [errortype char('0'+vars.warning)];
label_x = 'Number of Integration Points';
if isempty(vars.quantities{1})
    qs = fieldnames(RP.derivedQuantities); 
    qs = qs(ismember(qs,fieldnames(RP.contributions)));
else
    qs = vars.quantities;
end
if isempty(qs), fig = []; return; end
if ~vars.fignumber, vars.fignumber=2; end
fig = figure(vars.fignumber);
clf(vars.fignumber); fig.UserData = [];
ax = axes(fig); ax.NextPlot = 'add';
ax.TickLabelInterpreter = vars.opts{2}; ax.FontSize = vars.opts{4};
ax.YScale = 'log'; ax.Box = 'on';
xlabel(label_x,vars.opts{:})
ylabel(label_y,vars.opts{:})
colors = ax.ColorOrder;
markers = ['o';'s';'d';'^';'v';'>';'<'];
MC = {'Marker', 'Color'}; mc = {markers colors};
if length(qs)==1, MC = MC(end:-1:1); mc = mc(end:-1:1); end
lines = cell(1,length(qs)+length(vars.index)); idx = 1;
lines = lines(1:end-length(qs)==1-length(vars.index)==1);
qnames = cell(1,length(qs)); conv_ = NaN(1,0);
for it1 = 1:length(qs)
    n1 = MC{2}; v1 = mc{2}(idx,:);
    q = qs{it1}; 
    [conv,x] = RP.computeError(q,vars.index,errortype); 
    if isempty(conv), continue; end
    useBarPlot = (isscalar(vars.index)&&length(x{1})<3)||length([x{:}])<3;
    if useBarPlot, conv = [conv{:}]; end
    qnames{it1} = q;
    if length(qs)>1 && ~useBarPlot
        lines{idx} = plot(ax,NaN,NaN,'o-',n1,v1,'DisplayName',qnames{it1});
        idx = idx+1;
    end
    if it1==1
        xt = unique([x{:}]);
        xticks(xt(max(end-4,1):end));
    end
    for it2 = 1:length(vars.index)
        if it2>length(conv), break; end
        n2 = MC{1}; v2 = mc{1}(mod(it2-1,7)+1,:);
        if isempty(conv_) && useBarPlot
            conv_ = NaN(2,length(qs));
            names = cell(1,length(qs));
            x_ = [x{:}]; if isscalar(x_), x_ = [NaN x_(1)]; end
            conv_(end-1+isscalar(conv):end,1) = conv;
            names(1) = qnames(it1); ndx = 2;
            if length(vars.index)>1, break; end
        elseif useBarPlot
            conv_(end-1+isscalar(conv):end,ndx) = conv;
            names(ndx) = qnames(it1); ndx = ndx+1;
        elseif length(x{it2})==1 && it2==length(vars.index)
            if isempty(conv{end}), continue; end
            if isempty(conv_)
                conv_ = zeros(size(qnames)); x_ = x{it2}; 
                names = cell(size(qnames)); ndx = 1;
            end
            conv_(ndx) = conv{end}; names{ndx} = qnames{it1}; ndx = ndx+1; 
        else
            plot(ax,x{it2},conv{it2},'-',n1,v1,n2,v2);
        end
        if it1==length(qs) && length(vars.index)>1 && isempty(conv)
            if it1>1, v3 = 'k'; a = 'o'; else, v3 = v1; a = 'o-'; end
            if vars.index(it2)<length(RP.contours)
                nm = sprintf('Contour %d',it2);
            elseif vars.index(it2)==length(RP.contours)
                nm = 'Background Contour';
            elseif vars.index(it2)==length(RP.contours)+1
                nm = 'Expansion'; a = 'o';
            end
            lines{idx} = plot(ax,NaN,NaN,a,n1,v3,n2,v2,'DisplayName',nm);
            idx = idx+1;
        end
    end
end
if ~isempty(conv_)
    if ~isempty(ax.Children)
        subplot(1,2,1,ax);
        ax_ = subplot(1,2,2,'Parent',fig);
        ax.YLabel.String = [ax.YLabel.String ' (Integral)'];
        fig.Position(3) = 2.5*fig.Position(4);
    else
        ax_ = ax;
    end
    conv_(:,ndx:end) = []; names(ndx:end) = []; % remove empty content
    bar(ax_,x_,conv_)
    if any(isnan(x_)), ax_.XTick = x_(2); end
    ax_.FontSize = vars.opts{4};
    ax_.YScale = 'log';
    ylabel(ax_,label_y,vars.opts{:})
    ax_.YLim = [min(conv_,[],'all')*2e-1 max(conv_,[],'all')*12];
    legend(ax_,names,vars.opts{:}); legend(ax_,'boxoff')
    xlabel(ax_,label_x,vars.opts{:})
    ax_.TickLabelInterpreter = vars.opts{2};
    if length(vars.index)==2 && useBarPlot
        xtl = {'Integral' 'Expansion'};
        ax_.XTickLabel = xtl(~isnan(x_));
        ax_.XLabel.String = '';
    end
    if length(fig.Children)>2
        ax_.YLabel.String = [ax_.YLabel.String ' (Expansion)'];
    end
end
if ~isempty(lines)
    legend([lines{:}],'location',vars.legendLocation,vars.opts{:})
    legend(ax,'boxoff')
end
if any(vars.ylim), ylim(ax,vars.ylim); end
if any(vars.xlim), xlim(ax,vars.xlim); end
end
