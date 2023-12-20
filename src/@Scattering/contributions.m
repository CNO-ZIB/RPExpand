function out = contributions(sc, w_n, q)
%CONTRIBUTIONS evaluate residues based on QNMs (only for point source)

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

assert(sc.pointSource,'Only implemented for a single point source.');
assert(~isempty(sc.evsDir), 'You must provide resonance results.');
position = sc.keys.position; strength = sc.keys.strength;

nm = 0.5;
[~,loc] = ismember(w_n,sc.eigenvalues.'); % get location of eigenvalues
sc.quantities.(q).keys.fieldIds = loc; % specify ids for post processes

field = [sc.evsDir filesep 'fieldbag.jcm'];
ndx = [1:6 8:2:(2*length(w_n)+6)];
data = zeros(ndx(end),1); 
data(ndx) = [position(:);strength(:);w_n(:)];
if size(sc.resonance,1)~=ndx(end) || any(data(ndx)~=sc.resonance(ndx))
    E = sc.export({field},position);
    if length(loc)~=length(w_n), error('Eigenvalues do not match'), end
    E = nm*E(:,:,loc);
    E = reshape(sum(E.*strength(:),1),1,size(E,2),[]);
    w_n = reshape(w_n,1,1,[]);
    data(8:2:end) = E;
    sc.resonance = data;
    tags = {};
else
    E = sc.resonance(8:2:end);
    tags = sc.tags;
end

if strcmp(q,'DipoleEmission')
    out = -0.5.*1i*E.*E;
elseif sc.quantities.(q).logical(1)
    if isempty(tags), tags = sc.f({conj(w_n)}); sc.tags = tags; end
    tags = {{field},tags{1}};
    if length(tags{2})>length(tags{1})
        tags{1} = repmat(tags{1},1,length(tags{2}));
    end
    v = sc.f([tags;tags([2 1])],q); % run post processes
    nm = reshape(nm,1,[]);
    sz = [1 size(v{1})];
    v = {nm.*reshape(v{1},sz) conj(nm).*reshape(v{2},sz)};
    out = cat(3,1i*E.*v{1},conj(1i*E).*v{2});
else
    tags = {repmat({field},1,length(w_n))};
    v = sc.f(tags,q);
    out = reshape(nm,1,[]).*reshape(v{1},[1 size(v{1})]);
end
out = -2i*pi*out;
end
