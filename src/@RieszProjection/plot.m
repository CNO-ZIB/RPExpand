function varargout = plot(RP, varargin)
%PLOT plot results
%   You can plot custom quantities with the default plotting function if it
%   is scalar or provide your own plotting routine. You just need to set a 
%   function handle to the field 'plot' of the specified quantity. All
%   parameters passed to RieszProjection/plot will be passed as a struct.
%   If additional parameters are required you can edit the definition of
%   the parser. The expansion of the radiation pattern, provided with the
%   interface 'Scattering' to JCMsuite, has a custom plotting function.
%
%   PLOT(RP) A list of quantities is displayed which are available for
%   plotting. The left click with the mouse will plot the quantity under
%   the cursor.
%
%   PLOT(RP,quantity) Plots the provided quantity
%
%   PLOT(RP,quantity,index) Plots the data corresponding to index, which
%   can be an integer or an array of integers. The indices correspond to
%   eigenmodes and are displayed on the upper axis of the plot showing the
%   integration contour. If there are n modal contributions, the index n+1
%   gives the contribution of the background, n+2 the sum of the plotted
%   contributions and n+3 the reference solution. If you want the sum to
%   include all contributions consider the use of 'add2background' to
%   reduce the number of displayed modes. For error plots not the number of
%   eigenmodes but the number of contours is relevant. If there are n
%   contours, n+1 is the index that refers to a comparison with the
%   reference soltution.
%
%   PLOT(...,PARAM1,VAL1,PARAM2,VAL2,...) Specifies
%   parameters with name value pairs which include:
%      w0 (numeric): Sets the expansion points.
%      fignumber (integer): Sets the number of the figure.
%      ylim (numeric): Sets the ylim-property of the axis. If [0 0]
%         (default) the limits will be determined automatically.
%      xlim (numeric): Sets the xlim-property of the axis. If [0 0]
%         (default) the limits will be determined automatically.
%      interpreter (string): Defines the interpreter for all the text
%         displayed in the figures. The default is 'latex'.
%      fontSize (numeric): Defines the fontSize for all text displayed
%         in the figures. The default is 12.
%      legendLocation (string): Sets the location of the legend. The
%         default is 'northeast'.
%      xunit (string): When plotting an expansion the unit on the x-axis
%         can be changed from angular frequency (1/s, default) to a 
%         wavelength possible inputs are 'nm', 'um', 'mm' and 'm'.
%      yunit (string): When plotting an expansion a unit can be added to
%         the axis label. E.g. J in the case of ElectricFieldEnergy. It
%         will appear automatically in square brackets.
%      contourNumbers (logical): If true (default) on the upper axis of the
%         plot showing the integration contours, the indices corresponding
%         to the contours are displayed.
%      quantities (cell): Select which quantities you want to include to
%         the error plot.
%      add2bg (int): Integer accounting for modes which are not of
%         interest and can be added to the background.
%      errortype (string): Defines the error type. Possible values are
%         'relative' (default), 'absolute' and 'minimum'. The last option
%         displays the minimum of absolute and relative errror. 
%      custom (cell): May contain parameters for custom plot functions
%      data (cell): The expansion corresponding to w0.
%
%   fig = PLOT(...) Returns the figure handle.
%
%   The parameters you provide are saved such that you do not have to
%   provide them again if you replot the same quantity. Some are changed
%   for all quantities and others only for those they are relavant for. If
%   you do not provide a quantity and/or an index, these parameters will
%   be taken from the current figure. E.g. if you have plotted the
%   normalized decay rate with an instance 'rp' of the class 
%   'RieszProjection', you can change the unit of the x-axis from 1/s to nm
%   by simply typing: rp.plot('xunit','nm').
%
%   see also figure, legend, ylim, xlim

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: April-2023

% construct the parser if not yet there and parse the input arguments
persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/plot';
    addOptional(parser,'quantity','?',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addOptional(parser,'index',0,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'}));
    addParameter(parser,'w0',[],...
        @(x)validateattributes(x,{'numeric'},{}));
    addParameter(parser,'fignumber',0,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'}));
    addParameter(parser,'ylim',[0 0],...
        @(x)validateattributes(x,{'numeric'},...
        {'2d' 'numel' 2 'nondecreasing' 'real'}));
    addParameter(parser,'xlim',[0 0],...
        @(x)validateattributes(x,{'numeric'},...
        {'2d' 'numel' 2 'nondecreasing' 'real'}));
    addParameter(parser,'interpreter','latex',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'fontSize',12,...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'})); 
    addParameter(parser,'legendLocation','northeast',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'xunit','1/s',...
        @(x)validateattributes(x,{'char'},{}));
    addParameter(parser,'yunit','',...
        @(x)validateattributes(x,{'char'},{}));
    addParameter(parser,'contourNumbers',true,...
        @(x)validateattributes(x,{'logical'},{'scalar'})); 
    addParameter(parser,'quantities',{[]},...
        @(x)validateattributes(x,{'cell'},{'vector'})); 
    addParameter(parser,'add2bg',[],...
        @(x)validateattributes(x,{'numeric'},{'integer' 'nonnegative'}));
    addParameter(parser,'errortype','relative',...
        @(x)validateattributes(x,{'char'},{'nonempty'}));
    addParameter(parser,'custom',{},...
        @(x)validateattributes(x,{'cell'},{'vector'}));
    addParameter(parser,'data',[],...
        @(x)validateattributes(x,{'numeric'},{}));
end

% check which quantities can be plotted
qs = fieldnames(RP.derivedQuantities).'; 
qs_ = {'ComplexPlane'}; if ~isempty(qs), qs_ = ['Error' qs_]; end
ch = getChildren(RP.quantities,qs);
if ~isempty(ch), qs = [qs ch]; end
q_ref = fieldnames(RP.reference)';
if isempty(q_ref), q_ref = {}; end
idx = ismember(q_ref,[qs 'fields' 'frequencies']);
qs = [q_ref(~idx) qs];

% check which optional arguments are given
if nargin>1
    qs = [qs_ RP.availableQuantities];
    try quantity = validatestring(varargin{1},qs); 
    catch 
        quantity = []; 
        if ischar(varargin{1})
            try validatestring(varargin{1},parser.Parameters);
            catch
                error('Unknown quantity %s',varargin{1})
            end
        end
    end
    try fn = RP.plotSettings.fignumber.(quantity); catch, fn = 0; end
    fig = get(groot,'CurrentFigure'); % get current figure if exists
    if isempty(quantity) % add to current figure
        try qidx = fig.UserData; catch, qidx = {'?' 0}; end
        if isnumeric(varargin{1})
            varargin = [qidx(1) varargin];
        else
            varargin = [qidx varargin];
        end
    elseif ~isempty(fig) && fn~=fig.Number
    elseif nargin>2 && ~isnumeric(varargin{2})
        try idx = fig.UserData{2}; catch, idx = 0; end
        varargin = [quantity idx varargin(2:end)];
    else
        varargin{1} = quantity;
    end
end

% save name of the instance for convenient hrefs
if isempty(RP.name_)
    RP.name_ = inputname(1); 
end

parse(parser,varargin{:});
vars = parser.Results; 
q = vars.quantity; 
inside(RP,vars.w0) % expansion frequencies must be contained in the contour

vars.maxIndex = RP.nModes+3;
if strcmp(q,'ComplexPlane')
    vars.maxIndex = length(RP.contours); % passed to plotContours
    if RP.includeConjugatedPoles % RP.nModes is probably not yet defined
        vars.maxIndex = ceil(vars.maxIndex/2);
    end
elseif strcmp(q,'Error')
    vars.maxIndex = length(RP.contours)+1;
end

% check wether the index is valid
validateIndex(vars.index,vars.maxIndex,q)

% load and save settings for the plot
vars = loadAndSavePlotSettings(RP,q,vars,parser);
vars.opts = {'Interpreter' vars.interpreter 'FontSize' vars.fontSize};

% plot
switch lower(q)
    case 'complexplane'
        fig = plotContours(RP,vars);
    case lower(RP.availableQuantities)
        if ~vars.fignumber
            vars.fignumber = find(strcmp(RP.availableQuantities,q))+2;
        end
        if isfield(RP.quantities.(q),'plot')
            fig = RP.quantities.(q).plot(RP,vars);
        else, fig = plotQuantity(RP,vars);
        end
    case 'error'
        fig = plotError(RP,vars);
    otherwise
        allQs = [qs_ qs(cellfun(@(x)~RP.quantities.(x).hidden,qs))];
        href = '<a href="matlab: %s.plot(''%%s'')">%%s</a>\n';
        href = sprintf(href,RP.name_);
        str = ['You can plot the following '...
            'quantities:\n' repmat(href,1,length(allQs))];
        qs = cell(1,2*length(allQs));
        qs(1:2:end-1) = allQs; qs(2:2:end) = allQs;
        str = sprintf(str,qs{:});
        if ~strcmp(vars.quantity,'?')
            str = sprintf(['Unknown quantity %s. ' str], vars.quantity); 
        end
        str = str(1:end-1);
        disp(str);
        return;
end
if ishandle(fig), fig = {fig}; end
for it = 1:length(fig)
    fig{it}.UserData = {vars.quantity vars.index 'w0' vars.w0};
end
if ~isempty(fig), RP.plotSettings.fignumber.(q) = fig{1}.Number; end 
if nargout>0, varargout = fig; end

    function vars = loadAndSavePlotSettings(RP,q,vars,parser)
        usingDefaults = parser.UsingDefaults;
        parameters = parser.Parameters;
        pS = RP.plotSettings;
        if vars.index, pS.add2bg = []; end
        locals = {'fignumber','yunit','ylim', ...
                  'index','xlim','legendLocation'};
        ce = any(strcmp(q,{'ComplexPlane' 'Error'}));
        if ~any(strcmp(usingDefaults,'xunit'))
            try pS.xlim=rmfield(pS.xlim,'expansion'); catch, end
        end
        for it1 = 1:length(parameters)
            parameter = parameters{it1}; 
            if any(strcmp(parameter,{'data','w0'})), continue; end
            islocal = strcmp(parameter,locals);
            if ismember(parameter,usingDefaults)
                if strcmp(parameter,'legendLocation')
                    if strcmp(q,'ComplexPlane') % default for contours
                        vars.legendLocation = 'northwest';
                    end
                end
                if ~isfield(pS,parameter) % check if contained in pS
                    continue;
                end
                if ~any(islocal)
                    vars.(parameter) = pS.(parameter);
                elseif isfield(pS.(parameter),q)
                    vars.(parameter) = pS.(parameter).(q);
                elseif isfield(pS.(parameter),'expansion') && ~ce
                    vars.(parameter) = pS.(parameter).expansion;
                end
            else
                if strcmp(parameter,'quantity'), continue; end
                if ~any(islocal)
                    pS.(parameter) = vars.(parameter);
                elseif any(islocal(1:3)) || (any(islocal(4:6)) && ce)
                    pS.(parameter).(q) = vars.(parameter);
                elseif any(islocal(4:6))
                    pS.(parameter).expansion = vars.(parameter);
                end
            end
        end
        RP.plotSettings = pS;
        vars.warning = any(vars.index);
    end
end

function validateIndex(index,maxIndex,quantity)
%VALIDATEINDEX Check if the index is valid
%   As the contours are defined first, the number of contours does not
%   necessarily follow from the number of modes.

if isscalar(index) && ~index
elseif any(~index)
    error('All indices must be >0.')
elseif strcmp(quantity,'ComplexPlane') && any(index>maxIndex)
    error(['The index must be >0 and <=%d, which is '...
        'the index of the background contour'], maxIndex)
elseif strcmp(quantity,'Error') && any(index>maxIndex)
    error(['The index must be >0 and <=%d, which is the index '...
        'of the physical quantity'],maxIndex)
elseif any(index>maxIndex)
        error(['The index must be >0 and <=%d, which is '...
            'the index of the reference solution'],maxIndex)
end
end

function ch = getChildren(quantities,qs)
ch = cell(size(qs));
for it = 1:length(qs)
    if isfield(quantities.(qs{it}),'children')
        ch{it} = quantities.(qs{it}).children; 
    end
end
ch = [ch{:}];
if ~isempty(ch), ch = [ch getChildren(quantities,ch)]; end
end