function varargout = inside(RP,ps,shift)
%INSIDE Check if points are contained within the contour
%   in = INSIDE(RP,ps) Returns an array of logicals indicating whether the
%   corresponding entry of the input argument ps is contained in the
%   contour.
%
%   in = INSIDE(RP,ps,shift) If a shift is specified, the points must be
%   contained in the contour whose boundary is shifted in normal direction.
%
%   INSIDE(...) An error is thrown if any of the points is not contained in
%   the contour.
assert(isnumeric(ps))
if isempty(ps) || isempty(RP.contours)
    if nargout, varargout{1} = []; end
    return; 
end
if ~nargout
    if isreal(ps), ps = ps([1 end]);
    else
        ch = convhull(real(ps),imag(ps),'simplify',true);
        ps = ps(ch(1:end-1));
    end
end
sh = RP.shape{2}; if nargin==2, shift = [0 0]; end
if sh=='c' % circle
    cr = RP.shape{1};
    in = abs(ps-cr(1))/(cr(2)+shift(1))<=1;
elseif sh=='p' % polygon
    ps1 = RP.shape{4}; s = diff(ps1);
    n = s./abs(s)*exp(-1i*pi/2);
    a = abs(angle([n(2:end) n(1)].*exp(-1i*angle(n))));
    v = ps1(2:end); a_ = angle(n)+a/2;
    ps1 = v+1./cos(a/2).*(shift(1)*cos(a_)+1i*shift(2)*sin(a_));
    in = inpolygon(real(ps),imag(ps),real(ps1),imag(ps1));
elseif sh=='e' % ellipse
    A = RP.shape{4}; ps = (ps-A(1)-1i*A(2))*exp(-1i*A(5));
    a = real(ps)/(A(3)+A(4)+shift(1));
    b = imag(ps)/(A(3)-A(4)+shift(1));
    in = a.^2+b.^2<=1;
end
if nargout, varargout{1} = in; return; end
if ~all(in)
    error('All expansion points must be contained in the outer contour')
end
end