function [ps,rs,ds,p] = zerosPoles(RP,poles,varargin)
% Compute zeros, poles and their derivatives. This method does some pre
% and post processing. Please refer to the method 'locatePoles'.

persistent parser; 
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'RieszProjection/zerosPoles';
    addParameter(parser,'quantity','',...
        @(x)validateattributes(x,{'char'},{'vector'}));  
    addParameter(parser,'nMax',[],...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','integer'}));  
    addParameter(parser,'cutOff',[],...
        @(x)validateattributes(x,{'numeric'},{'nonnegative','scalar'}));
    addParameter(parser,'plot',[],...
        @(x)validateattributes(x,{'logical'},{'scalar'}));
end

% default cutoff and nMax
args = struct('cutOff',RP.cutOff,'nMax',RP.nMax);
% default return values
ps = []; rs = []; ds = []; p = parser;

if isempty(RP.contours)
    error('Contour required!')
elseif isempty(varargin)
    printPossibilities(RP,poles); return; 
else 
    args.quantity = validatestring(varargin{1},RP.availableQuantities); 
end

for it = 2:length(varargin)
    if isfloat(varargin{it}) && varargin{it}<1
        args.cutOff = varargin{it};
    elseif isfloat(varargin{it}) && varargin{it}<Inf
        args.nMax = round(double(varargin{it}));
    elseif ischar(varargin{it})
        for it_ = it:2:length(varargin)
            args.(varargin{it}) = varargin{it+1};
        end
    end
end
parse(parser,args); q = parser.Results.quantity;
cutoff = parser.Results.cutOff; nmax = parser.Results.nMax;
assert(cutoff<1)
while isfield(RP.quantities.(q),'parents')
    q = RP.quantities.(q).parents{1};
end
if ~isfield(RP.derivedQuantities,q)
    RP.derivedQuantities.(q) = {}; RP.up2date(1) = false;
end

if ~RP.up2date(1), RP.evaluate; end
if poles && isempty(RP.defineBg)
    ps = cell(1,2);
    if ~isempty(RP.selectedPoles)
        ps{1} = reshape(RP.poles([RP.selectedPoles{:}]),[],1);
    end
    if ~isempty(RP.expansionPoints)
        ps{2} = RP.expansionPoints([1 end]); 
    end
    RP.defineBg = unique(cat(1,ps{:}));
end

args = {cutoff nmax q};
if poles && RP.nPoints
    results = cell(3,RP.nModes*(1+RP.quantities.(q).quadratic));
    for it = 1:size(results,2)
        ndx = mod(it-1,RP.nModes)+1;
        args{2} = max(2,length(RP.selectedPoles{ndx}));
        [results{:,it}] = RP.locatePoles(it,args{:});
    end
    ps = [results{1,:}]; ds = [results{3,:}]; rs = [results{2,:}];
else
    % negative index for zeros and positive index for poles
    [ps,rs,ds] = RP.locatePoles((-1+2*poles)*length(RP.contours),args{:});
end
end


function printPossibilities(RP,poles)
qs = RP.availableQuantities;
qs(cellfun(@(x)RP.quantities.(x).hidden,qs)) = [];
for it = 1:length(qs)
    while isfield(RP.quantities.(qs{it}),'parents')
        qs{it} = RP.quantities.(qs{it}).parents{1};
    end
end
qs = unique(qs);

f = 'computePoles'; if ~poles, f = 'computeZeros'; end
href = sprintf('<a href="matlab: %s.%s(''%%s'')">%%s</a>\n',RP.name_,f);
str = ['You can use the following '...
    'quantities:\n' repmat(href,1,length(qs))];
q = cell(1,2*length(qs));
q(1:2:end-1) = qs; q(2:2:end) = qs;
str = sprintf(str,q{:});
str = str(1:end-1);
disp(str);
end
