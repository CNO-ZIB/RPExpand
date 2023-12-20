function addQuantity(sc,name,varargin)
%ADDQUANTITY add quantity based on given post process
%   With the method 'addPostProcess' you add a quantity that is marked
%   inactive and will not be passed to the class 'RieszProjection' via the
%   method 'getQuantities' unless it is explicitely added as a quantity.
%   This can be done using the name of the post process. In this case, all
%   quantities that will defined based on the same post process will
%   inherit its properties unless they are explicitly set. Of course, it is
%   also possibilbe to assign a quantity with a given set of keys and
%   parameters a new name. Further quantities can than be either based on
%   the orginal post process with default properties or the already defined
%   quantity.
%
%   ADDQUANTITY(sc,name) convert the post process to a quantity of the same
%   name with default keys. 
% 
%   ADDQUANTITY(sc,name,keys) change the keys for parameter
%   substitution associated with the quantity with the specified name
%
%   ADDQUANTITY(sc,name,parent,keys) copy the quantity with the
%   specified name assign it new keys for parameter substitution and save
%   it using the specified alias
%
%   ADDQUANTITY(...,PARAM1,VAL1,PARAM2,VAL2,...) add fields to the
%   specified quantity. These may be used to evaluate the target quantity
%   from the results of the post-process. For all quantities the following
%   properties can be set: 
%      quadratic (logical): If true the quantity is regarded quadratic in
%         the electric field.
%      hidden (logical): If true the quantity is marked as a hidden
%         quantity.
%      evaluate (function_handle): Defines a custom function used to
%         extract the quantity of interest from the result of the post
%         process. It must take the name of the quantity and tags, which 
%         refer to the resultbag, and return the corresponding values of 
%         the target quantity in a double array of size (d,n), where d is
%         a custom data dimension and n the number of tags.
%      m (double): derivatives and vector-valued quantities must share the
%         same axis. For this reason, m gives the index to the last data
%         dimension, i.e., v(m+1:end) are partial derivatives of the data
%         points.
%      reference (cell): If for the reference solution another post
%         process is used, keys and parameters can be provided in a cell
%         array. It must have the form {keys,PARAM1,VAL1,PARAM2,VAL2,...}
%         or {PARAM1,VAL1,PARAM2,VAL2,...}.
%
%   see also quantity, addPostProcess, RieszProjection, getQuantities

assert(isempty(sc.parent.f),sprintf(['After a call of the method '...
    '''getQuantities'', which is called in the constructor\nof the '...
    'class ''RieszProjections'',the quantities must not be changed.']))

% check if parent or keys are provided
pk = {name struct}; m = 1; n = 0;
n = n + (~isempty(varargin) && ...
    (ischar(varargin{1}) && ...
    ismember(varargin{1},fieldnames(sc.quantities))));
pk(1:n) = varargin(m:n); m = m+n;
n = n + (~isempty(varargin(m:end)) && isstruct(varargin{m}));
pk(2:1+(m==n)) = varargin(m:n);
[parent,ks] = pk{:};

if mod(length(varargin(n+1:end)),2)
    error('Please check the arguments'); 
end

q = sc.quantities.(parent);
q.logical(2:3) = [false true]; % customized quantities are activated
sc.quantities.(name) = setParameters(q,ks,varargin(n+1:end));
end

function q = setParameters(q,ks,ps)
reference = 0;
for it = 1:2:length(ps)
    switch lower(ps{it})
        case 'quadratic'
            q.logical(1) = ps{it+1};
        case 'hidden'
            q.logical(2) = ps{it+1};
        case 'reference'
            reference = it;
        otherwise
            q.(ps{it}) = ps{it+1};
    end
end
for f = fieldnames(ks).'
    q.keys.(f{1}) = ks.(f{1});
    if isfield(q,'reference') && ~reference
        q.reference.keys.(f{1}) = ks.(f{1});
    end
end
q = editQuantity(q);
if reference
    if isfield(q,'reference'), ref = q.reference; else, ref = q; end
    ks = struct; ps = ps{reference+1};
    if isstruct(ps{1}), ks = ps{1}; ps = ps(2:end); end
    q.reference = setParameters(ref,ks,ps);
elseif isfield(q,'reference')
    q.reference = editQuantity(q.reference);
end
end
