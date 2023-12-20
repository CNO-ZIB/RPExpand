function varargout = computePoles(RP,varargin)
%COMPUTEPOLES compute the poles (eigenvalues)
%   zs = COMPUTEPOLES(RP,...) If no additional arguments are provided, a 
%   list of quantities that can be used to compute the poles of the system
%   is displayed. Otherwise, all the poles of the specified quantity and 
%   inside the contour are located. The defaults of the additional
%   arguments are are defined by properties of the class RieszProjection 
%   with the same name.
%   Arguments:
%     quantity (char) - compute the poles of the specified quantity
%     cutOff (double, <1) - poles with normalized residues smaller
%     than the specified cutoff are discarded (a smaller number will
%     result in a larger number of poles)
%     nMax (integer, >=1) defines the size of the Hankel matrix and
%     should be larger then the expected number of eigenmodes
%
%   see also: RieszProjection/cutOff, RieszProjection/nMax

% save name of the instance for convenient hrefs
if isempty(RP.name_), RP.name_ = inputname(1); end
[ps,rs,ds,parser] = zerosPoles(RP,true,varargin{:});
if isempty(varargin), return; end
if nargout, varargout{1} = ps; end

q = parser.Results.quantity;
if ~isempty(ds)
    RP.derivatives.poles = ds(1);
    RP.derivatives.residues.(q) = ds(2);
end

ndx = inside(RP,ps);
if RP.quantities.(q).quadratic
    ndx_ = imag(ps)<=0; ps = ps(ndx_); ndx = ndx_ & ndx;
else, ps = sort(ps,'ComparisonMethod','real');
end
RP.poles_ = ps;

RP.selectedPoles_ = num2cell(find(ismember(RP.poles_,ps(ndx))));
RP.nModes = length(RP.selectedPoles_);

RP.up2date(2:3) = [true false];
RP.contributions = struct;
if RP.quantities.(q).m == 1
    RP.contributions.(q).residues = reshape(rs,[1 size(rs)]);
end

% update plots
RP.up2date(5) = false;
updatePlots(RP,[RP.availableQuantities 'Error' 'ComplexPlane'])
end

