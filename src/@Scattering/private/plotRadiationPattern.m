% plot the radiation pattern either 2d or 3d

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022
function fig = plotRadiationPattern(RP,vars)

% get the data and the keys to reconstruct the mesh in theta and phi 
rad = []; ref = [];
hz = strcmpi(vars.xunit,'1/s'); % if x-unit is not 1/s revert order
q = RP.f.quantities.(vars.quantity).keys;
if ~vars.index, index = 1:RP.nModes+3; else, index = vars.index; end
if isempty(vars.w0), vars.w0 = RP.referencePoints; end

if any(index<RP.nModes+3)
    rad = vars.data; args = {vars.w0 'precision', Inf};
    if isempty(rad), rad = RP.computeExpansion(vars.quantity, args{:}); end
    if ~hz, rad = rad(end:-1:1,:,:); end
end
if any(index==RP.nModes+3)
    ref = RP.computeReference(vars.quantity,vars.w0);
    if ~hz, ref = ref(end:-1:1,:); end
end
theta = q.gridPointsTheta; phi = q.gridPointsPhi;

index(ismember(index,vars.add2bg(vars.add2bg~=RP.nModes+1))) = []; 

if length(phi)==2
    [theta,idx] = sort([-theta theta]);
    rad = rad(:,idx,:); ref = ref(:,idx);
    [theta,idx] = unique(theta); 
    rad = rad(:,idx,:); ref = ref(:,idx);
end

% get normalization
tot = sum(rad,3); if isempty(ref), ref = tot; end
m = size(tot,1); idx = zeros(1,RP.nModes+1,'logical');
idx(end) = true; idx(vars.add2bg) = true; bg = sum(rad(:,:,idx),3);
nm_ = max(ref(:)); rad = cat(3,rad(:,:,1:end-1),bg,tot,ref); 

% plot
fig = figure(vars.fignumber); clf; xL = [Inf -Inf]; yL = [Inf -Inf];
is3D = all([length(theta) length(phi)]>2); ax = cell(1,m); 
if is3D
    [Phi,Theta] = meshgrid(deg2rad(phi),pi/2-deg2rad(theta));
    sz = size(Theta); zL = [Inf -Inf];
else
    theta = deg2rad(theta)+pi/2; % 2d
end
ax_ = subplot(1,m+1,1); disableDefaultInteractivity(ax_) % for legend
ax_.CLim = [1 RP.nModes+3];
for it1 = 1:m
    ax{it1} = subplot(1,m+1,it1+1); ax{it1}.DataAspectRatio = [1 1 1];
    ax{it1}.NextPlot = 'add'; ax{it1}.Visible = 'off';
    % in case that the axis is made visible
    ax{it1}.TickLabelInterpreter = vars.interpreter;
    ax{it1}.FontSize = vars.fontSize; ax{it1}.CLim = [1 RP.nModes+3];
    ax{it1}.Clipping = 'off'; ax{it1}.Projection = 'perspective';
    % plot
    
    for idx = index
        if is3D
            R = real(reshape(rad(it1,:,idx)./nm_,sz));
            [X,Y,Z] = sph2cart(Phi,Theta,R);
            s = mesh(ax{it1},X,Y,Z); s.CData = idx*ones(sz);
            s.FaceAlpha = 0.1; s.FaceColor = 'flat';
        else
            [X,Y] = pol2cart(theta,real(rad(it1,:,idx))./nm_);
            l = plot(ax{it1},X,Y);
        end
        if idx==RP.nModes+2 && is3D
            s.FaceColor = 'none';
        elseif idx==RP.nModes+3
            if is3D
                s.EdgeColor = 'none'; s.FaceAlpha = 0.5;
            else
                l.LineStyle = ':'; l.LineWidth = 2;
            end
        end
        xL = [min(xL(1),min(X(:))) max(xL(2),max(X(:)))];
        yL = [min(yL(1),min(Y(:))) max(yL(2),max(Y(:)))];
        if is3D, zL = [min(zL(1),min(Z(:))),max(zL(2),max(Z(:)))]; end
    end
end

% Link the aspect ratios, reduce distances and put labels

ax = [ax{:}]; 
if is3D
    link = linkprop(ax,{'CameraPosition','CameraUpVector',...
        'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(fig, 'StoreTheLink', link);
    ax(end).ZLim = zL; ax(end).YLim = yL; ax(end).XLim = xL;
    ax(end).CameraPosition = [2.5 -2.5 2.5];
    ax(end).CameraTarget = [0 0 zL(2)/2];
    ax(end).CameraUpVector = [0 0 1];
else
    linkaxes(ax);
    ax(end).YLim = yL; ax(end).XLim = xL;
end
opts = ['Units','normalized','Visible','on',vars.opts];
[x,xu] = unit(RP.referencePoints,vars.xunit);
x = sort(x);
if strcmp(xu,'nm'), pt = '%.0f'; else, pt = '%.3g'; end
sh = 1/(10*m); % move to correct positions
for it1 = 1:m
    pos = it1/m-1/m/2-0.03*is3D;
    text(ax_,pos,0.1,sprintf(pt,x(it1)),opts{:});
    ax(it1).Position(1:3) = [sh+(it1-1)*(1/m-sh/m) 0.15 1/m-2*sh/m];
end
if hz, wl = '$\omega_0$ '; else, wl = '$\lambda_0$ '; end
if strcmpi(vars.interpreter,'latex')
    t = [wl '[$\mathrm{' xu '}$]'];
else
    wl = split(wl,{'$'}); t = [wl{2} ' [' xu ']'];
end
text(ax_,0.47,0.05,t,opts{:});

ax_.Visible = 'off'; ax_.NextPlot = 'add';
ax_.Position = [sh-0.02 0 1-2*(sh-0.02) 0.95];
ax_.XTick = 0.5:1:m; ax_.XLim = [0 m]; ax_.XTickLabels = string(x);

% legend
leg = cell(1,length(index));
for it1 = 1:length(index)
    idx = index(it1);
    if is3D
        s = mesh(ax_,zeros(2),zeros(2),zeros(2)); 
        s.CData = ones(2)*idx; s.FaceAlpha = 0.1;
        s.FaceColor = 'flat';
    else
        l = plot(ax_,NaN);
    end
    if idx<RP.nModes+1
        leg{it1} = sprintf('Mode %d',idx);
    elseif idx==RP.nModes+1
        leg{it1} = 'Background';
    elseif idx==RP.nModes+2
        leg{it1} = 'Sum';
        if is3D, s.FaceColor = 'none'; end
    elseif idx==RP.nModes+3
        leg{it1} = 'Reference';
        if is3D
            s.EdgeColor = 'none'; s.FaceColor = 'flat';
            s.FaceAlpha = 0.5; s.CData = s.CData+1;
        else
            l.LineWidth = 2; l.LineStyle = ':';
        end
    end
end
legend(ax_,leg,vars.opts{:},'location',vars.legendLocation);
theta = q.gridPointsTheta;
split = ~isempty(vars.custom) && vars.custom{1};
if split && ~is3D&&(max(theta)<=90 || min(theta)>=90)
    upper = max(theta)<=90;
    figs = split2DPattern(fig,upper); 
    fig = [{fig},figs];
end
end

function figs = split2DPattern(fig,upper)
% sT is the number of segments in theta and should be a multiple of four
% sR is the number of segments in r direction
% RLim is the maximal radius and sT*N the number of points for a circle
sT = 12; sR = 2; RLim = 1.2; N = 20;
figs = cell(1,length(fig.Children)-2);
circs = (RLim/sR:RLim/sR:RLim).*exp((2i*pi*(0:sT*N-1).'/(sT*N)));
S.Vertices = [[real(circs(:)) imag(circs(:))];0 0]; 
S.Faces = NaN(sT*sR,2*N+2);
S.FaceVertexCData = zeros(sT*sR,1);
for it = 1:sT*sR
    if it<=sT
        k = it*N+1-(it*N==sT*N)*sT*N;
        S.Faces(it,1:N+2) = [(it-1)*N+1:it*N k size(S.Vertices,1)];
    else
        s = S.Faces(it-sT,1:N+1);
        S.Faces(it,:) = [s+sT*N s(end:-1:1)];
    end
    if upper
        S.FaceVertexCData(it) = mod(it-1,sT)<sT/2;
    else
        S.FaceVertexCData(it) = mod(it-1,sT)>=sT/2;
    end
end
fig1 = figure(123); clf; fig1.Visible = 'off'; ax_c = axes(fig1); 
ax_c.DataAspectRatioMode = 'manual'; ax_c.Visible = 'off'; 
ax_c.NextPlot = 'add'; ax_c.Colormap = [0.8 0.8 0.8; 1 1 1];
p = patch(ax_c,S,'FaceColor','flat'); 
p.LineWidth = 0.5; p.EdgeAlpha = 1; ax_c.CLim = [0 1];
p.EdgeColor = [0.15 0.15 0.15]; p.FaceColor = 'flat'; p.FaceAlpha = 1;
 
args = {'Interpreter','latex','FontSize',fig.Children(1).FontSize};
args1 = [args 'VerticalAlignment' 'cap'];
txt = text(0,RLim,'$0^o$','HorizontalAlignment','center',args{:});
h = txt.Extent(end); txt.Position(2) = txt.Position(2)+h*0.4;
text(0,-RLim-0.07*h,'$0^o$','HorizontalAlignment','center',args1{:});

for it = 1:sT/4
    theta = it*2*pi/sT; t = it*360/sT; r = RLim;
    lable = sprintf('$%d^o$',t); lable_ = ['$-' lable(2:end)];
    args1 = {'VerticalAlignment' 'bottom' 'HorizontalAlignment' 'left'};
    t = text(r*sin(theta),r*cos(theta),lable,args{:},args1{:});
    args1{end} = 'right'; r = RLim - 0.2*h;
    t1 = text(-r*sin(theta),r*cos(theta),lable_,args{:},args1{:});
    if it==sT/4
        t.String = ['$\pm' t.String(2:end)];
        t1.String = ['$\mp' t1.String(3:end)];
        t.VerticalAlignment = 'middle'; t1.VerticalAlignment = 'middle';
        t.HorizontalAlignment = 'left'; t1.HorizontalAlignment = 'right';
        t.Position = t.Position + h*0.2*[sin(theta) cos(theta) 0];
        t1.Position = t1.Position - h*0.4*[sin(theta) cos(theta) 0];
    else
        args1{end-2} = 'top'; args1{end} = 'left'; r = RLim-0.1*h;
        text(r*sin(theta),-r*cos(theta),lable_,args{:},args1{:});
        args1{end-2} = 'top'; args1{end} = 'right'; r = RLim;
        text(-r*sin(theta),-r*cos(theta),lable,args{:},args1{:});
    end
end

args2 = [args 'VerticalAlignment' 'middle','HorizontalAlignment','center'];
for it = 1:length(fig.Children)-2
    ax = fig.Children(it); scaling = 0;
    for line = ax.Children(end:-1:1).'
        s = round(max(sqrt(line.XData.^2+line.YData.^2)),2);
        if s<scaling, continue; else, scaling = 100*s; end
        rest = mod(scaling,sR);
        if rest~=0, scaling = scaling+sR-rest; end
        scaling = scaling/100;
    end
    fig_ = figure(fig.Number+it); clf(fig_);
    ax_ = copyobj(ax_c,fig_);
    theta = pi/sT;
    for r = RLim/sR:RLim/sR:RLim
        if upper, sgn_u = 1; else, sgn_u = -1; end
        lable = sprintf('$%.2f$',r*scaling);
        text(ax_,sgn_u*r*sin(theta),sgn_u*r*cos(theta),lable,args2{:});
        lable = sprintf('$-%.2f$',r*scaling);
        text(ax_,-sgn_u*r*sin(theta),-sgn_u*r*cos(theta),lable,args2{:});
    end
    copyobj(ax.Children,ax_)
    for line = ax_.Children(1:length(ax.Children)).'
        line.XData = line.XData/scaling; line.YData = line.YData/scaling;
    end
    ax_.FontSize=11;
    leg = legend(ax_.Children(length(ax.Children):-1:1),args{:});
    leg.String = fig.Children(end-1).String; leg.NumColumns = 2;
    if upper
        leg.Location = 'southeast';
    else
        leg.Location = 'northeast';
    end
    figs{it} = fig_;
end
close(fig1);
end
