function varargout = locatePoles(RP,index,cutoff,nmax,q)
%LOCATEPOLES get poles based on values along a closed contour

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: April-2023

n = abs(index); m = size(RP.derivedQuantities.(q){n},2);
dvs = m-RP.quantities.(q).m; q2 = RP.quantities.(q).quadratic;

center = sum(RP.contours{n}(:))/numel(RP.contours{n});
normx = round(max(abs(RP.contours{n}(:)-center)),4,'significant');
nEvs = min(nmax*(1+q2),floor(numel(RP.contours{n})/2));
mu = zeros(2*nEvs,m); 
w = reshape(RP.contours{n},[1 1 size(RP.contours{n})]);
vars = struct('quantity',q,'dvs',dvs);
for it1 = 1:2*nEvs
    g = ((w-center)./normx).^(it1-1);
    mu(it1,:) = RP.quad(g,index,vars);
end

if ~dvs, mu = sum(mu,2); end
normy = mu(1);
center = center/normx; mu = mu./normy;

H = hankel(mu(1:nEvs,1),mu(nEvs:end-1,1));
H2 = hankel(mu(2:nEvs+1,1),mu(nEvs+1:end,1));

[R,D] = eig(H2,H,'vector');
V = (D.^(0:nEvs-1)).';
rs = diag(R.'*H*R)./sum(V.'.*R.',2).^2;
selection = abs(rs./max(rs))>cutoff; 
[ps,ndx] = sort(normx*(D(selection)+center).','ComparisonMethod','real');
rs = rs(selection)*normy; rs = rs(ndx).';

% match the number of poles with positive and negative imaginary parts 
% if H is a real matrix, its eigenvalues must be symmetric with respect
% to the real axis
if max(abs(imag(H)))/max(abs(H)) < 1e-5 && ~RP.nPoints
    n = abs(sum(imag(ps)>0)-sum(imag(ps)<0));
    [~,ndx] = sort(abs(rs./max(rs)));
    rs(ndx(1:n)) = []; ps(ndx(1:n)) = [];
end
if q2
    ndx = imag(ps)<=0; ps = [ps(ndx) ps(~ndx)]; rs = [rs(ndx) rs(~ndx)];
end
ndx = inside(RP,ps); ps_ = ps(ndx); rs = rs(ndx);

% solve for derivatives if given
derivatives = struct([]);
    
if dvs
    nDerivatives = size(mu,2); 
    if ismember('derivativeNames',properties(RP.f))
        derivativeNames = RP.f.derivativeNames;
    else
        derivativeNames = split(strtrim(sprintf('dp%d ',1:nDerivatives)));
    end
    rsn = rs/normy; n = length(ps_);
    nk = (-1:2*n-1).'; wnk = (ps_/normx-center).^nk; mu = mu(1:2*n,:);
    nm = regexp(derivativeNames,'(?<=_)\w+','match','once');
    for it = 2:nDerivatives
        if n==1
            d = (mu(2,it)-mu(2,1)*mu(1,it));
        elseif n>1
            d = [nk(2:end).*wnk(1:end-1,:).*rsn wnk(2:end,:)]\mu(:,it);
        else
            d = [];
        end
        derivatives(1).(nm{it}) = d(1:n)*normx;
        derivatives(2).(nm{it}) = d(n+1:end)*normy;
    end
end

results = {ps rs derivatives};
varargout(1:nargout) = results(1:nargout);
end