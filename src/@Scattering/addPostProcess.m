function addPostProcess(sc,name,varargin)
%QUANTITY add quantity based on a post process
%   If the quantity already exists it will be replaced. Some default
%   quantities are added automatically during construction (please refer 
%   to the function private/addDefaultQuantities). For these default
%   quantities some parameters are set before they are added. If you want
%   to replace one of these quantities please check if you need to edit the
%   function private/editQuantity.
%   
%   QUANTITY(sc,name) The quantity of the given name based on a post
%   process ['postprocesses' filesep name '.jcmpt'].
%
%   QUANTITY(sc,name,...) you can provide the following additional 
%   arguments in an arbitrary order: 
%      keys (struct): keys for parameter substitution in the input file of
%         the post process, passed to jcmwave_solve
%      postprocess (char): file for the postprocess ending with .jcmpt that
%         can either be a file name located in the directory postprocesses
%         or the absolute file name.
%      reference (cell): a cell array containing arguments, as described in
%         this section, for the reference solution. This can be needed, if
%         the reference solution is based on different post process or 
%         requires special keys.
%      logical (logical): The argument 'logical' provides answers to the
%         binary questions:
%         logical(1): is it a quadratic quantity? (default false)
%         logical(2): is it a hidden quantity? (default false)
%         logical(3): is this quantity active? (default false)
%         Only quantities with logical(3) = true are passed to the
%         RieszProjection class with the method 'getQuantities'. The
%         logical may be scalar or only have two elements. In this case,
%         the remaining entries are set to their defaults. 
%      evaluate (function_handle): Sometimes the quantity of interest is 
%         not directly returned by the post process and some further post
%         processing is required. A corresponding file handle can be
%         provided, that takes the name of the quantity and tags, which 
%         refer to the resultbag, and returns the target quantity as double 
%         array of size (d,n), where d is a custom data dimension and n the
%         number of tags.
%      m (double): derivatives and vector-valued quantities must share the
%         same axis. For this reason, m gives the index to the last data
%         dimension, i.e., v(m+1:end) are partial derivatives of the data
%         points.
%
%   see also customizeQuantity, Scattering

% delete old definition and set some defaults
if isfield(sc.quantities,name)
    sc.quantities.(name) = rmfield(sc.quantities,name); 
end
sc.quantities.(name) = struct('keys',struct,'logical',false(1,3), ...
                              'evaluate',@extractResult,'m',1);

% set quantity
sc.quantities.(name) = setPP(sc.quantities.(name),name,varargin);

    function out = setPP(out,name,args)
        postProcess = name; reference = 0;
        out.id = name; out.parent = sc; 
        for it = 1:length(args)
            switch class(args{it})
                case 'struct'
                    out.keys = args{it};
                case 'char'
                    postProcess = args{it};
                case 'cell'
                    reference = it;
                case 'logical'
                    if isscalar(args{it}), args{it}(2) = false; end
                    out.logical(1:length(args{it})) = args{it};
                case 'function_handle'
                    out.evaluate = args{it};
                case 'double'
                    out.m = args{it};
                otherwise
                    error('addQuantity: unexpected argument')
            end
        end
        if reference
            args = args{reference}; ref = rmfield(out,'logical');
            out.reference = setPP(ref,['reference' name],args);
        end
        out.jcmpt = checkPP(postProcess);
    end

    function pP = checkPP(pP)
        exclude = {'DipoleEmission' 'Superposition' 'PointEvaluation'};
        if ismember(pP,exclude), pP = []; return; end
        if exist(pP,'file'), return; end
        pP = [sc.postproDir filesep pP];
        if ~strcmp(pP(end-6:end),'.jcmpt')
            pP = [pP '.jcmpt'];
        end
        if ~exist(pP,'file')
            error('Scattering:addQuantity:missingPostProcess',...
                'Please provide the file %s',pP)
        end
    end
end
