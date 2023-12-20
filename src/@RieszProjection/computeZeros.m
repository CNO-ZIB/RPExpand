function varargout = computeZeros(RP,varargin)
%COMPUTEZEROS compute the zeros of the specified quantity
%   zs = COMPUTEZEROS(RP,...) If no additional arguments are provided, a
%   list of quantities whose zeros can be evaluated is displayed.
%   Otherwise, all the zeros of the specified quantity and inside the
%   contour are located. The defaults of the additional arguments are
%   defined by properties of the class RieszProjection with the same name.
%   Arguments:
%     quantity (char) - compute the zeros of the specified quantity
%     cutOff (double, <1) - zeros with normalized residues smaller
%     than the specified cutoff are discarded (a smaller number
%     will result in a larger number of zeros)
%     nMax (integer, >=1) defines the size of the Hankel matrix and
%     should be larger then the expected number of zeros
%
%   see also: RieszProjection/cutOoff, RieszProjection/nMax

% save name of the instance for convenient hrefs
if isempty(RP.name_), RP.name_ = inputname(1); end
[zs,~,ds,parser] = zerosPoles(RP,false,varargin{:});
if isempty(varargin), return; end
plt = parser.Results.plot; q = parser.Results.quantity;

RP.zeros.(q) = zs; if nargout, varargout{1} = zs; end
if ~isempty(ds), RP.derivatives.zeros.(q) = ds(1); end
plt = (isempty(plt) && ~nargout) || (~isempty(plt) && plt); 

if isempty(zs) || ~plt, return; end

try fig = figure(RP.plotSettings.fignumber.ComplexPlane);
catch
    fig = RP.plot('ComplexPlane');
end
ax = fig.Children(end);
if contains(ax.Children(1).DisplayName,'z')
    delete(ax.Children(1));
end
dName = '\omega_z';
if strcmpi(fig.Children(1).Interpreter,'latex')
    dName = ['$' dName '$'];
end
ax.NextPlot = 'add';
plot(ax,zs,'bo','DisplayName',dName)
updatePlots(RP,[RP.availableQuantities 'Error']);
end

