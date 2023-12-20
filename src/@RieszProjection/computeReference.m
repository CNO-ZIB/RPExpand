function varargout = computeReference(RP,quantity,w0)
%COMPUTEREFERENCE Compute the reference solutions on the real line
%   COMPUTEREFERENCE(RP, quantity) Computes the reference solution of
%   the given quantity.
%   The results are saved in the property 'reference', which is a struct
%   with field names corresponding to the quantities.
%
%   COMPUTEREFERENCE(RP, quantity, w0) The reference solutions are
%   computed at the given real frequencies.
%
%   reference = COMPUTEREFERENCE(...) The reference solution is
%   returned. Its shape is (m,n) with m being the number of reference
%   points and n the data dimension. If derivatives are available
%   reference(:,1) may be the scalar quantity and reference(:,2:end) the
%   derivatives with respect to different design parameters. 

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022

% make sure w0 is given
if nargin<3, w0 = RP.referencePoints; end
if isempty(w0), return; end
q = RP.quantities.(quantity);
if isfield(q,'parents')
    qs = cell(size(q.parents));
    for it = 1:length(qs)
        qs{it} = RP.computeReference(q.parents{it},w0);
    end
    if nargin(q.evaluate)>length(qs), qs{end+1} = w0; end
    reference = q.evaluate(qs{:});
else
    reference = referenceSolution(RP,quantity,w0);
    if isfield(q,'evaluate')
        reference = {reference w0};
        reference = q.evaluate(reference{1:nargin(q.evaluate)});
    end
end
if nargout, varargout{1} = reference; end

    function reference = referenceSolution(RP,quantity,w0)
        w0 = [real(w0(:)) imag(w0(:))];
        w = RP.reference.frequencies;
        w1 = uniquetol([w;w0],'ByRows',true);
        
        % if necessary evaluate the scattering problems at the expansion
        % points, query only those which are not already present
        [s,x] = ismembertol(w1,w,'ByRows',true);
        if isempty(RP.reference.fields)
            fs = RP.f({w0(:,1)+1i*w0(:,2)}); fs = fs{1};
            RP.reference.fields = fs;
            RP.reference.frequencies = w0;
        elseif sum(s) ~= length(w1)
            fs = cell(1,size(s,1));
            fs(s) = RP.reference.fields(x(s));
            fs_ = RP.f({w1(~s,1)+1i*w1(~s,2)}); fs(~s) = fs_{1};
            RP.reference.frequencies = w1;
            RP.reference.fields = fs;
            
            % update frequency indices of existing quantities
            [~,x] = ismembertol(w,w1,'ByRows',true);
            names = fieldnames(RP.reference);
            for qn = names(3:end).'
                ndx = zeros(1,size(w1,1),'logical');
                ndx(x(RP.reference.(qn{1}).w0)) = true;
                RP.reference.(qn{1}).w0 = ndx;
            end
        else
            fs = RP.reference.fields;
        end
        
        s = ismembertol(w1,w0,'ByRows',true);
        
        % check if existing results are up to date and can be reused
        w0_q_old = false(1,size(w1,1)); values_old = [];
        if isfield(RP.reference,quantity)
            w0_q_old = RP.reference.(quantity).w0;
            values_old = RP.reference.(quantity).values;
        end
        w0_q_new = w0_q_old | s.';
        up2date = w0_q_old(w0_q_new);
        if any(~up2date)
            if any(up2date)
                values(up2date,1:size(values_old,2)) = values_old;
            end
            qn = ['reference' quantity];
            values_new = RP.f({fs(w0_q_new&(~w0_q_old))},qn);
            values(~up2date,1:size(values_new{1},1)) = values_new{1}.';
            RP.reference.(quantity).values = values;
            RP.reference.(quantity).w0 = w0_q_new;
        else
            values = RP.reference.(quantity).values;
        end
        reference = values(s(w0_q_new),:);
    end
end