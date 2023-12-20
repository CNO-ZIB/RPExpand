function viewEigenmodes(sc,p,c)
%VIEWEIGENMODES show the eigenmodes, if available, in JCMview
%   If the property 'evsDir' has been set, the fieldbag it contains is
%   opened. Otherwise, if the parent is set and eigenmodes are available,
%   the residues corresponding to them are computed with a contour integral
%   and for each eigenmode the corresponding fieldbag is opened. 
%
%   VIEWEIGENMODES(sc) open all available eigenmodes with JCMview.
%
%   VIEWEIGENMODES(sc,p) show the eigenmodes associated with the indices p
%   in the contour with index one of the corresponding instance of the 
%   class 'RieszProjection'. If nPoints>0, i.e., if small contours exist,
%   the first eigenvalues of the contours with indices p are shown. 
%
%   VIEWEIGENMODES(sc,p,c) show the eigenmodes associated with the indices
%   p{n} in the contour with index c(n) of the corresponding instance of 
%   the class 'RieszProjection'. In this case p must be a cell array with
%   length(c) elements. 

if ~isempty(sc.evsDir)
    jcmwave_view([sc.evsDir filesep 'fieldbag.jcm']);
end

if isempty(sc.parent.fields), fprintf(1,'No eigenmodes available'); end
rp = sc.parent;

if nargin<3
    if nargin==1, p = 1:rp.nModes; end; c = 1;
    selection = rp.selectedPoles;
    if rp.nPoints
        c = p; p = cellfun(@(x){1:length(x)},rp.selectedPoles);
    else
        selection = {[selection{:}]};
    end
end
if ~iscell(p), p = {p}; end
res = cell(2,sum(cellfun(@length,p))); n = 1;
for it1 = 1:length(c)
    tags = rp.fields{c(it1)}; ndx = p{c(it1)};
    for it2 = ndx
        weights = reshape(rp.getWeights(it2,c(it1)),1,[]);
        wn = rp.poles(selection{c(it1)}(it2));
        res(:,n) = sc.sum(tags,weights,wn); n = n+1;
    end
end
job_ids = [res{1,:}]; modes = res(2,:);
if ~isempty(job_ids), [~] = jcmwave_daemon_wait(job_ids); end
for mode = modes
    jcmwave_view(mode{1});
end

