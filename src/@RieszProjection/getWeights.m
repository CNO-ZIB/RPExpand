function weights = getWeights(RP,p,c,w0)
%GETWEIGHTS get weights for quadrature
%   The contour integral I is computed with the quadrature
%   I = sum(weights.*values,[3 4]) or r_i = sum(weights.*values,3)
%   for Gauss-Kronrod and the rapezoidal rule, respectively.
%
%   weights = GETWEIGHTS(RP,p) the weights for the residue of the pole with
%   index p inside the contour with index one are returned. 
%   If p = 0, all the poles are considered, i.e., the integral I is the sum
%   of the contributions of all the poles inside the contour.
%
%   weights = GETWEIGHTS(RP,p,c) the weights for the residue of the pole
%   with index p inside the contour with index c are returned. 
%
%   weights = GETWEIGHTS(RP,p,c,w0) the facotr 1/(w-w0) is multiplied to
%   the weights.
%
%   The contributions to a quantity are 1/(2i*pi)*(I_1 + I_2 + ... + I_bg)
%   where the indices of the integrals refer to the different poles. 
%
%   see also: RieszProjection

if nargin<4, w0 = []; end; if nargin<3, c = 1; end 
vars = struct('quantity','weights','p',p);
if ~p
    g_p = 1;
elseif c==length(RP.fields)
    g_p = RP.poles([RP.selectedPoles{:}]); % consider all poles
else
    g_p = RP.poles(RP.selectedPoles{c});
end
weights = RP.quad(g_p,c,vars);
if ~isempty(w0) % 1/(w-w0) is multiplied 
    w = reshape(RP.contours{n},[1 1 size(RP.contours{n})]);
    weights = weights./(w-w0(:)); 
end
end

