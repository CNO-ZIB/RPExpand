function varargout = computeExpansion(RP, varargin)
%COMPUTEEXPANSION Perform a modal expansion
%   COMPUTEEXPANSION(RP) All quantities available for an expansion are
%   listed and can be selected with the mouse. A left click will call this
%   method with default parameters or parameters that had been used when
%   the selected quantity has been expanded previously. 
%
%   COMPUTEEXPANSION(RP, quantity) An expansion of the given quantity is
%   evaluated. The interface that has been passed to the constructor must
%   be able to derive this quantity from the solutions of the linear
%   systems. The argument 'quantity' can either be a string or a cell array 
%   containing strings. If the property 'precision' is set, the number of
%   points is refined adaptively until the maximal error is sufficiently
%   small or the number of points on the corresponding contour is larger
%   than the value of the keyword argument 'upperbound', whose default is
%   1024. The property 'precision' can be either a scalar or a vector with
%   two elements. In the latter case, the first value is used as the bound 
%   of the absolute error and the second one as the bound of the relative
%   error. For absolute and relative error control you can set the first or 
%   second value to zero, respectively. If you pass the keyword argument 
%   precision, it is used instead of the property. If the precision is Inf
%   no refinements will be made.
%   The residues of the quantity are evaluated and can be accessed with
%   rp.contributions.(quantity).residues. Furthermore, the expansion is
%   plotted. For scalar quantities a plotting function is available. If the
%   quantity has a field 'plot' containing a function handle, the custom
%   function is used instead. 
%
%   COMPUTEEXPANSION(RP,quantity,w0) Use the expansion frequencies w0
%   instead of the property 'expansionPoints'.
%   
%   COMPUTEEXPANSION(...,PARAM1,VAL1,PARAM2,VAL2,...) The above arguments
%   can be followed by key-value pairs, which include:
%      precision (numeric): Determines when to stop the adaptive
%         refinement. The default is the value of the property of the same
%         name, which is Inf by default.
%      upperbound (numeric): Fixes the maximal number of integration points
%         a single contour can have. The default is 1024.
%      plot (logical or cell): Plot the expansion with default arguments
%         (logical) or with the arguments defined in the cell array. If
%         nargout==0 the default is true and false otherwise. 
%
%   expansion = COMPUTEEXPANSION(...) the expansion is returned which is
%   of size (k,m,n) with k being the number of expansion points, m is the
%   data dimension (usually one) and n is the number of considered
%   eigenmodes  plus one. The physical quantity is therefore
%   q = sum(expansion,3);. 
%
%   see also RieszProjection RieszProjection/plot

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: April-2023

persistent parser;
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/computeExpansion';
    parser.StructExpand = false;
    addOptional(parser,'quantity','?',@validatequantity);
    addOptional(parser,'w0',[],...
        @(x)validateattributes(x,{'numeric'},{'vector'}));
    addParameter(parser,'precision',0,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'}));  
    addParameter(parser,'upperbound',2^10,...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','integer'}));
    addParameter(parser,'plot',false,...
        @(x)validateattributes(x,{'logical','cell'},{'vector'}));
end

parse(parser,varargin{:});
vars = parser.Results;
vars.usingDefaults = parser.UsingDefaults;
q = vars.quantity; if ischar(q), q = {q}; end

% check if arguments for plot are provided
args = {}; plt = true;
if iscell(vars.plot), args = vars.plot; else, plt = vars.plot; end
plt = plt || ~nargout;

% save name of the instance for convenient hrefs
if isempty(RP.name_), RP.name_ = inputname(1); end

% Check if any expansion points are defined
w0 = {vars.w0(:) RP.expansionPoints RP.referencePoints};
w0 = [w0(~cellfun(@isempty,w0)) {[]}]; vars.w0 = w0{1};
if isempty(vars.w0), error('No expansion points defined'); end

% If no quantity is specified print links to available quantities
if q{1}(1)=='?', printPossibilities(RP), return; end

varargout = cell(1,nargout); qs = RP.availableQuantities;
for it = 1:length(q)
    quantity = validatestring(q{it},qs,mfilename,'quantity');
    % Get modal contributions and background
    expansion = get_expansion(RP,quantity,vars);
    if plt, RP.plot(q{it},args{:},'w0',vars.w0,'data',expansion); end
    if it<=nargout, varargout{it} = expansion; end
end
if plt, qs = qs(~ismember(qs,q)); end
updatePlots(RP,[qs 'Error']);
end

% expand single quantity
function out = get_expansion(RP,qn,vars)
q = RP.quantities.(qn);
if isfield(q,'parents')
    nParents = length(q.parents);
    out = cell(1,nParents);
    for it = 1:nParents
        out{it} = get_expansion(RP,q.parents{it},vars);
    end
    if nargin(q.evaluate)>nParents, out{end+1} = vars.w0; end
    out = q.evaluate(out{:});
    return
end
args = {'precision',vars.precision,'upperbound',vars.upperbound};
if ~isfield(RP.contributions,qn) 
    if ~isfield(RP.derivedQuantities,qn)
        RP.derivedQuantities.(qn) = {}; RP.up2date(1) = false;
    end
end
if q.quadratic==2 && ismember('w0',vars.usingDefaults) ...
        && ~isempty(RP.referencePoints)
    vars.w0 = RP.referencePoints;
end
bg = RP.computeContributions(qn,vars.w0,args{:});
% Get reference solution if reference points are defined
RP.computeReference(qn)
% the second axis is the data axis
residues = RP.contributions.(qn).residues;
residues = reshape(residues,1,size(residues,2),[],1+q.quadratic);
wn = reshape(RP.poles([RP.selectedPoles{:}]),1,1,[]);
if q.quadratic, wn = cat(4,wn,conj(wn)); end
out = cat(3,sum(-1/(2i*pi)*residues./(wn-vars.w0),4),bg);

% combine modal contributions if desired
if size(out,3)>RP.nModes+1
    it1 = 0;
    for it2 = 1:length(RP.selectedPoles)
        it1 = it1+1;
        if length(RP.selectedPoles{it2})==1 && it1~=it2
            out(:,:,it2) = out(:,:,it1);
        else
            s = it1; it1 = s-1+length(RP.selectedPoles{it2});
            out(:,:,it2) = sum(out(:,:,s:it1),3);
        end
    end
    out = out(1:it2);
end
if isfield(q,'evaluate')
    out = {out vars.w0}; 
    out = q.evaluate(out{1:nargin(q.evaluate)}); 
end
end

% Print all quantities that can be expanded to the command window
function printPossibilities(RP)
qs = RP.availableQuantities;
qs(cellfun(@(x)RP.quantities.(x).hidden,qs)) = [];
href = '<a href="matlab: %s.computeExpansion(''%%s'')">%%s</a>\n';
href = sprintf(href,RP.name_);
str = ['You can expand the following '...
    'quantities:\n' repmat(href,1,length(qs))];
q = cell(1,2*length(qs));
q(1:2:end-1) = qs; q(2:2:end) = qs;
str = sprintf(str,q{:});
str = str(1:end-1);
disp(str);
end

function validatequantity(quantity)
if ischar(quantity) || iscellstr(quantity) %#ok<ISCLSTR>
    return;
end
error(['The input parameter quantity is expected to be a string' ...
    'or a cell array containing strings.'])
end
