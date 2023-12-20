function [err,x] = computeError(RP, quantity, index, kind)
%COMPUTEERROR compute an estimate of the integration error
%   [err,x] = COMPUTEERROR(RP, quantity, index, kind) The errors of the
%   integrals along contours corresponding to the specified indices are
%   returned. The kind of error can be 'absolute' ('abs'), 'relative'
%   ('rel'), or 'minimum' ('min') which is the minimum of absolute and
%   relative error weighted with absolute and relative error tolerances
%   defined in the property 'precision'.
%   The returned values are cell arrays of shape (1,m) if the argument
%   'index' has m elements. If a Gauss-Kronrod quadrature rule is used for
%   integration, the error of the background contour is given by the
%   difference between Gauss and Gauss-Kronrod quadrature.
%   The indices can take values between 1 and n+1 where n is the number of
%   contours. The error of the full expansion with respect to the reference
%   solution is err{k} if index(k) = n+1. The output x{k} contains in this 
%   case the total number of integration points (from all contours). The
%   output err{k} is a vector and err{k}(m) is the error of the integration
%   with x{k}(m) abscissae.
%   The default kind of error is minimum and the default index is n+1.
%
%   err = COMPUTEERROR(...) only the error bound(s) of the best available
%   approximation(s) is returned. A scalar or a vector (if length(index>1)
%   is returned.

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

% check which data is available
dat = isfield(RP.derivedQuantities,quantity);
refD = isfield(RP.reference,quantity);
nContours = length(RP.contours); 
tol = RP.precision_; out_bool = false; % set defaults

% set defaults 
if nargin<4, kind = 'r'; if nargin<3, index = nContours+1; end; end
% define tolerance for absolute and relative errors, out => print warnings
if isnumeric(kind), tol = kind; kind = 'm'; out_bool = true; out = false;
elseif length(kind)>1, out = logical(kind(2)-'0'); kind = kind(1); 
end

if ~tol(1), kind = 'r'; elseif ~tol(end), kind = 'a'; end
if any(isinf(tol)) && kind(1)=='m', kind = 'r'; end
% check if there is enough data to compute an error
cnt = ~isempty(RP.contours); % check if any contours exist
gk = cnt && RP.shape{3}>1 && length(RP.weights{end})>1; % Gauss-Kronrod
cnt = cnt && (gk||~mod(size(RP.contours{end},1),2)); % error of contours
msg = [];
if ~(dat && (cnt||refD))
    msg = ['\nThere is not enough data for an error estimation of ' ...
        'the background contribution to the %s.\n'];
    index(nContours) = [];
elseif ~cnt && any(index~=nContours+1)
    msg = ['\nOnly the error with respect to the reference solution '...
        'is available for the quantity %s.\n'];
    index = nContours+1;
elseif ismember(nContours+1,index) && ~refD
    msg = ['\nThe implemented estimation of the error bound of the '...
        'physical quanity ''%s'' requires a reference solution. \n'...
        'Yet, estimations of error bounds of individual '...
        'contours are available.\n'];
    index(index==nContours+1) = [];
end

if ~isempty(msg) && out, fprintf(msg,quantity); end
if isempty(index), err = {}; x = {}; return; end

index(index==RP.nModes+1) = nContours;
index(index==RP.nModes+2) = nContours+1;

% if necessary update scattering solutions and derived quantities
if ~RP.up2date(1), RP.evaluate; end
% compute the error 
err = cell(1,length(index)); x = cell(1,length(index));
for it1 = 1:length(index)
    ndx = index(it1); gk_ = gk&&ndx==length(RP.contours);
    if ndx == nContours+1
        w0 = RP.referencePoints;
        ref = RP.computeReference(quantity);
        vs = RP.computeExpansion(quantity,w0,'precision',Inf);
        m = RP.quantities.(quantity).m;
        vs = vs(:,1:m,:); ref = ref(:,1:m);
        err{it1} = real(ref - sum(vs,3));
        x{it1} = sum(cellfun(@numel,RP.contours));
    else
        ref = RP.quad(1,ndx,struct('quantity',quantity));
        n = 1;
        if ~gk_ && nargout>1
            nv = numel(RP.contours{ndx})/2;
            while mod(nv,2)==0 && nv>2, nv = nv/2; n = n+1; end
            x{it1} = nv*2.^(n-1:-1:0);
        elseif nargout>1
            x{it1} = numel(RP.contours{ndx});
        end
        err{it1} = RP.quad(1,ndx,struct('quantity',quantity,'err',n));
        if gk_ && nargout>1, err{it1} = sum(err{it1},4); end
    end
    err{it1} = maxError(err{it1},ref);
end
if nargout==2, return; end
if all(cellfun(@isscalar,err)) || out_bool, err = [err{:}]; end

    function err = maxError(e_abs,ref)
        tol_ = tol;
        e_abs = abs(e_abs);
        e_rel = e_abs./abs(ref);
        if kind(1)=='m'
            [a,b] = min(cat(3,e_abs/tol(1),e_rel/tol(end)),[],3);
            if out_bool, tol_ = ones(size(tol)); end
            err = max(a.*reshape(tol_(b),size(b)),[],[1 2]);
        elseif kind(1)=='r'
            err = max(e_rel,[],[1 2]);
        elseif kind(1)=='a'
            err = max(e_abs,[],[1 2]);
        end
        err = squeeze(err).';
        if out_bool
            if length(err)>1
                v = RP.shape{1}; sp = RP.shape{2};
                if any(sp=='ec'), v = v(3:end); end
                h = abs(diff(v)); pathlen = sum(h);
                err = err < h/pathlen; 
            else
                err = (err<1);
            end
        end
    end
end
