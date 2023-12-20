function fig = viewFields(RP,vars)
%VIEWFIELDS visualisation of the total field
%   fig = VIEWFIELDS(RP,vars) Visualizes the expansion of the field. This
%   function is called by RieszProjection/plot and the following custom
%   parameters are available:
%      cmap (numeric): Colormap as returned by colormap(...). The colormap
%         is used for the visualization with matlab (see imshow). The
%         default is a red blue colormap.
%      component (char): defines the quantity to be displayed. The 
%         parameter 'component' can take the values real*, imag*, abs and 
%         log. The asterix stands for x, y and z. Be aware that the 
%         intensities will not sum up to the reference solution as taking 
%         the absolute value is not a holomorphic operation. If the
%         parameter component is not passed to the plot function, the
%         component that comes with the quantity is used. 
%      normalize (logical): If true, the values will be normalized to
%         range between -1 and 1. The default is true if the component is
%         neither abs nor log.
%      displayRange (numeric): Two-element vector [LOW HIGH] that controls
%         the display range. The default is [-1 1] if normalize = true
%         and [min max] otherwise.
%      d (numeric): A scalar that gives the width of a colored frame which
%         helps to identify the rows.
%
%DEPENDENCIES: 
%   This method requires the class RieszProjection
%
%   see also RieszProjection/plot imshow


persistent parser
if isempty(parser)
    red = [zeros(4,1);linspace(0,1,28).';ones(27,1);linspace(1,0.8,5).'];
    blue = [linspace(0.8,1,5).';ones(27,1);linspace(1,0,28).';zeros(4,1)];
    green = [zeros(4,1);linspace(0,1,28).';linspace(1,0,28).';zeros(4,1)];
    cm = [red green blue];
    parser = inputParser;
    parser.FunctionName = 'plotFields';
    addParameter(parser,'cmap',cm,...
        @(x)validateattributes(x,{'numeric' 'char'},{})); 
    addParameter(parser,'component','',...
        @(x)validateattributes(x,{'char'},{'vector'})); 
    addParameter(parser,'normalize',true,...
        @(x)validateattributes(x,{'logical'},{'scalar'})); 
    addParameter(parser,'displayRange',[0 0],...
        @(x)validateattributes(x,{'numeric'},...
        {'size', [1 2], 'increasing'}));
    addParameter(parser,'d',3,...
        @(x)validateattributes(x,{'numeric'},{'scalar','nonnegative'}));
end
quantity = RP.f.quantities.(vars.quantity); field = [];

parse(parser,vars.custom{:}); p = parser.Results;
if ismember('component',parser.UsingDefaults)
    p.component = quantity.component;
end

if ~vars.index, index = 1:RP.nModes+3; else, index = vars.index; end
if isempty(vars.w0), vars.w0 = RP.referencePoints; end
hz = strcmpi(vars.xunit,'1/s'); % if x-unit is not 1/s revert order

if any(index<RP.nModes+3)
    field = vars.data; args = {vars.w0 'precision', Inf};
    if isempty(field)
        field = RP.computeExpansion(vars.quantity, args{:}); 
    end
    sel = ismember(index,vars.add2bg); index(sel) = [];
    if ~isempty(vars.add2bg)
        sel = sum(field(:,:,vars.add2bg),3);
        field(:,:,RP.nModes+1) = field(:,:,RP.nModes+1)+sel;
    end
    field = field(:,:,index(index<RP.nModes+2));
end
tot = []; ref = [];
if any(index==RP.nModes+2), tot = sum(field,3); end
if any(index==RP.nModes+3)
    ref = RP.computeReference(vars.quantity,vars.w0);
end
field = cat(3,field,tot,ref);

wl = '\omega_0';
if ~hz, field = field(end:-1:1,:,:); wl = '\lambda_0'; end
sz = quantity.size(quantity.size>1);
field = reshape(field(:,:,end:-1:1),size(field,1),sz(1),sz(2),3,[]);
c = p.component(1); ndx = int16(p.component(end))-119;
switch c
    case 'r'
        field = real(field(:,:,:,ndx,:));
    case 'i'
        field = imag(field(:,:,:,ndx,:));
    case 'l'
        field = log(vecnorm(field,2,4));
    otherwise
        field = vecnorm(field,2,4);
end
field = permute(field,[3 5 2 1 4]);
sz = size(field,1:4);
field = reshape(field,sz(1)*sz(2),sz(3)*sz(4));

fig = figure(vars.fignumber);
clf(fig); ax = axes(fig,'Visible','off');
if ~any(p.displayRange), p.displayRange = []; end
if p.normalize && all(c~=['a' 'l'])
    m = min(field,[],'all'); M = max(field,[],'all');
    mM = max(-m,M);
    field = field/mM;
    if isempty(p.displayRange), p.displayRange = [-1 1]; end
end
args = {'Colormap',p.cmap,'DisplayRange',p.displayRange,'Parent',ax};
imshow(field, args{:}); ax.YDir = 'normal';
ax.NextPlot = 'add'; ax.DataAspectRatio = [1 1 1];
nFields = length(index); nFrequencies = length(vars.w0);
height = size(field,1); width = size(field,2); d = p.d;
x = {[0 width width 0] [1*d width-d width-1*d d]};
y = [[height height]-height/nFields [height height]];
y = {y [y(1:2)+d y(3:end)-d]};
name = cell(1,nFields); ndx = 1; line = cell(1,nFields);
args = {'FaceColor',p.cmap(floor(end/2),:),'LineWidth',1.5,...
    'FaceAlpha',1}; c = ax.ColorOrder;
for it = index
    if ~d, break; end
    if it<RP.nModes+1
        name{ndx} = sprintf('Mode %d',it);
    elseif it==RP.nModes+1
        name{ndx} = 'Background';
    elseif it==RP.nModes+2
        name{ndx} = 'Sum';
    elseif it==RP.nModes+3
        name{ndx} = 'Reference Solution';
    end
    p = polyshape(x,y); y = cellfun(@(x){x-height/nFields},y);
    ps = plot(ax,p,'FaceColor',c(mod(it-1,7)+1,:),'EdgeColor','none');
    line{ndx} = plot(ax,polyshape,args{:}); % for legend
    line{ndx}.EdgeColor = ps.FaceColor; ps.FaceAlpha = 1;
    ndx = ndx+1;
end
ax.XLim = [0,width]; ax.YLim = [0,height];
args = {'FontSize',vars.fontSize,'Interpreter',vars.interpreter};
if nFields>=nFrequencies
    args_ = {'NumColumns',1,'Location','WestOutside'};
else
    args_ = {'NumColumns',2,'Location','NorthOutside'};
end
legend([line{:}],name,args{:},args_{:});
if nFields>=nFrequencies, ax.OuterPosition = [0 0 0.9 0.9]; end
if isinf(vars.w0(1)), hold off, return; end
d = width/nFrequencies; xy = [d/2 0]; 
[x,xu] = unit(vars.w0,vars.xunit); x = sort(x);
if strcmp(xu,'Hz'), xu = ['P' xu]; x = x*1e-15; end
if strcmp(xu,'nm'), pt = '%.0f'; else, pt = '%.2g'; end
if strcmp(vars.interpreter,'latex')
    wl = ['$' wl '$']; xu = ['$' replace(xu,'s','\mathrm{s}') '$'];
end
for it = 1:nFrequencies
    t = text(xy(1),xy(2),sprintf(pt,x(it)),args{:});
    t.Position(2) = t.Position(2)-t.Extent(end); xy(1) = xy(1)+d;
    t.Position(1) = t.Position(1)-t.Extent(3)/2;
end
t_label = text(width/2,-t.Extent(end),[wl ' [' xu ']'],args{:});
t_label.Position(2) = t_label.Position(2)-t_label.Extent(end);
t_label.Position(1) = t_label.Position(1)-t_label.Extent(3)/2;
hold off
end
