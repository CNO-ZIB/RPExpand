% plot an expansion

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: July-2022
function fig = plotQuantity(RP, vars)
q = []; q_ref = []; fig = [];
requiredq = [true true];
if any(vars.index)
    requiredq(1) = any(vars.index<RP.nModes+3);
    requiredq(2) = any(vars.index==RP.nModes+3);
else
    vars.varning = false;
end

% get data
if isempty(vars.w0), vars.w0 = RP.expansionPoints; end
if requiredq(1)
    q = vars.data; args = {vars.w0 'precision' Inf};
    if isempty(q), q = RP.computeExpansion(vars.quantity, args{:}); end
end
if requiredq(2)
    w0_ref = [];
    if ~isempty(RP.referencePoints)
        w0_ref = RP.referencePoints;
    elseif vars.warning
        fprintf('No reference solution available.\n');
    end    
    if ~isempty(w0_ref)
        xref = unit(w0_ref,vars.xunit); 
        q_ref = RP.computeReference(vars.quantity,w0_ref);
        q_ref = q_ref(:,1);
    end
end

% set indices
emptyq = [isempty(q) isempty(q_ref)]; if all(emptyq), return; end
emptyq = emptyq & requiredq; M = RP.nModes+3;
if vars.index, inds = vars.index; else, inds = 1:M-(~requiredq(2)); end
if emptyq(1),inds=inds(inds==M); elseif emptyq(2), inds=inds(inds<M); end

% get ylabel from quantity name
a = find(ismember(vars.quantity,'A':'Z'));
a = [1 a(a>1) length(vars.quantity)+1]; 
a = a([true diff(a(1:end-1))>1 true]);
qname = join(arrayfun(@(x,y)vars.quantity(x:y),a(1:end-1),a(2:end)-1,...
    'UniformOutput',false));
% unitcomversions
[x,xu] = unit(vars.w0(:),vars.xunit);
if contains(xu,'s')
    wl = '$\omega_0$ '; 
else
    wl = '$\lambda_0$ ';
end
% plot
leg = cell(1,length(inds)-length(vars.add2bg(vars.add2bg~=RP.nModes+1)));
it2 = 1; fig = figure(vars.fignumber);
clf(fig)
ax1 = axes(fig); 
ax1.TickLabelInterpreter = vars.opts{2}; ax1.FontSize = vars.opts{4};
hold on;
if any(inds<M)
    final = zeros(size(q,1),1);
    for it1 = 1:RP.nModes
        if ismember(it1,vars.add2bg)
            q(:,1,end) = q(:,1,end)+q(:,1,it1); continue;
        elseif ~ismember(it1,inds)
            continue; 
        end
        final = final + real(q(:,1,it1));
        plot(x, real(q(:,1,it1)), '.-')
        leg{it2} = sprintf('Mode %d',it1); it2=it2+1;
    end
end
if ismember(RP.nModes+1,inds) && ~isnan(q(1,1,end))
    final = final + real(q(:,1,end)); 
    if length(x)==1, ln = '.'; else, ln = '--'; end
    plot(x, real(q(:,1,end)), ln); leg{it2} = 'Background'; it2=it2+1;
end
if ismember(RP.nModes+2,inds) && sum(~cellfun(@isempty,leg))>1
    plot(x, final, '.-'); leg{it2} = 'Sum';
end
if ismember(M,inds)
    plot(xref, real(q_ref), 'o');
    leg{end} = 'Reference';
end
hold off
if any(vars.ylim), ylim(vars.ylim); end
if any(vars.xlim), xlim(vars.xlim); end
if strcmpi(vars.interpreter,'latex')
    xlabel(ax1,[wl '[$\mathrm{' xu '}$]'],vars.opts{:})
    if ~isempty(vars.yunit)
        vars.yunit = ['$\mathrm{' vars.yunit '}$'];
    end
else
    wl = split(wl,{'$'});
    xlabel(ax1,[wl{2} ' [' xu ']'],vars.opts{:})
end
if ~isempty(vars.yunit), qname = [qname{1} ' [' vars.yunit ']']; end
ylabel(qname,vars.opts{:})
leg = leg(~cellfun(@isempty,leg));
if length(leg)<2, return; end
legend(leg,'location',vars.legendLocation,vars.opts{:}); legend('boxoff')
end
