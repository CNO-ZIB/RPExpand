%% Problem description
% We want to solve the Schroedinger equation in 1D for a quantum particle 
% of effective mass m* = 1 given the potential V(x) = -10 for -L<=x<=L and 
% 0 otherwise. Furthermore, h-bar is set to sqrt(2) and L = pi/sqrt(2). A
% thorough description of the physical problem and its discretization
% can be found in the paper 'An eigenvalue method for open-boundary quantum 
% transmission problems' by Shao et al. (1995, doi: 10.1063/1.360132). The
% numerical constants are taken from Gavin et al.: 'FEAST eigensolver for
% nonlinear eigenvalue problems' (2018, doi: 10.1016/j.jocs.2018.05.006).
% 
% From the physical point of view, we are interested in the question, of how
% the modes of the open system contribute to the transition amplitude at a
% potential well.
%
% This simple problem runs very fast and is well suited to get familiar
% with the capabilities of the class RieszProjection.
%
% This example covers two typical use cases. First, we assume that we are
% interested in a small energy interval, where we expect less than ten
% modes to dominate the spectrum. This allows for solving the eigenvalue 
% problem and performing the expansion with solutions of linear systems
% along a single contour. In the second case, the eigenvalues are given
% from an external source and the remaining task is the evaluation of the
% expansion. Here, the number of modes that are taken into account can be 
% much larger than in the first case. 

%% Prepare the workspace

% In most cases, it is nice to start with a clean workspace.
clearvars; clc

% In order to get reproducible results, the random number generation is
% reset before running this script.
rng(15)

% If your current folder is not the directory containing the class
% definition '@RieszProjection', you should add its parent directory to the
% Matlab path.
% addpath('/path/to/your/copy')
addpath('../src');

%% Construct an instance

% As we want to test the eigensolver, which is expected to work well for a
% small number of eigenvalues (about 5), we do not need to pass eigenvalues
% from an external source to the constructor. 
% The physical quantity we are interested in is the transmission which is 
% quadratic in the solutions of the Schroedinger equation. Therefore, the 
% conjugated poles must be included.
parameters.includeConjugatedPoles = true;

% Furthermore, define a range of frequencies (the expansion points), which
% represent the energy interval we are interested in. 
parameters.expansionPoints = 1:0.1:6;

% Define the number of integration points on the background contour
parameters.nPointsB = 128;

% and their minimal distance to the expansion points. In the case of an
% ellipse, the first value will be added to the mayor and the second value
% to the minor semiaxis. If the shape is selected to be a polygon, the first 
% and second values correspond to shifts in the x and y direction, respectively.
% If you want the outer contour to be a circle, you can provide a scalar
% value which is added to the radius.
parameters.minD = [0.4 2]; 
parameters.radius = 1e-2;

% As we expect to find poles below the real axis, we should add enough 
% space between expansion points and contour in the horizontal direction 
% and we must guess how many poles the contour will enclose.
parameters.nMax = 12; % Maximum number of eigenvalues

% Eigenvalues whose normalized residues are smaller than the cut off value
% are discarded.
parameters.cutOff = 1e-3; 

% Eventually, we define some reference points which allow us to check if
% the result converges to the expected solution. 
parameters.referencePoints = parameters.expansionPoints(1:2:end);

% As we provide a custom interface, we must tell which quantities can be
% expanded and which of them are linear. A linear quantity is required to
% solve the eigenvalue problem and if RieszProjection knows the available
% quantities, it will suggest them if you call the method 'expand' without 
% arguments.
qs.RightBoundary = struct('quadratic',false,'hidden',true);
qs.Transmission = struct('quadratic',true);
qs.Reflection = struct('quadratic',true);
parameters.quantities = qs;

% Now, an instance of the class RieszProjection can be constructed. Apart
% from the parameters we also hand a function handle to the constructor,
% which will be used to solve the linear systems, i.e., the Schroedinger
% equations, at the integration points and to derive the target quantities.
% It is defined at the bottom of this script.
rp = RieszProjection(@f, parameters);

%% Get the contours and expand

% First, we must construct the contours. We want an elliptic shape and the
% trapezoidal rule, which is the default.
rp.getContours

% Now, we can plot the contour and, subsequently, expand the target quantity.
% As our target quantity is a quadratic form, but the property 'poles' is
% empty, the linear quantity 'RightBoundary' is expanded automatically in 
% order to solve the eigenvalue problem. The solutions of the linear
% systems at the integration points are subsequently used to evaluate the
% expansions.
f1 = rp.plot('ComplexPlane');
rp.computeExpansion('Transmission');

% In order to see if the error decreases, we first plot the error. In
% contrast to other methods, we do not need the eigenvectors to estimate
% the error. As we are particularly interested in the expansion, we compare
% the result solutions of the linear system at the reference points, i.e.,
% solutions of the Schroedinger equation at the real line.
f2 = rp.plot('Error');

% The error is indeed decreasing and the expansion, based on the 
% eigenvalues, displayed in figure 1, matches the reference solution.
% There are two ways to further improve the accuracy: On the one hand, we 
% can simply increase the number of integration points. On the other hand, 
% since it is known that the convergence rate depends on the distance of 
% the contour to the closest pole, we could first define a new shape. 
% As the contour passes closely to one of the poles, this might indeed be 
% a good idea. We will do both and compare the results.

%% Refine the contour

% As we wish to reuse the existing solutions of the linear system, we
% double the number of integration points twice. Subsequently, we just
% expand the quantity we are interested in again. As the existing poles are
% used as initial guesses for the eigenvalues, it might matter if the number
% of integration points is increased at once or in two subsequent steps.
% For simplicity, we will just do the first. As the eigenvalues have not
% been set by the user, they will be recalculated automatically if expand
% is called after the property 'nPointsB' has been increased. If you want
% to improve the precision of the eigenvalues from an external source, you
% can call expand with the additional argument 'recalculateEvs'. E.g.,
% rp.expand('Transmission','recalculateEvs',true)
rp.nPointsB = 256;

% The figure containing the contours is updated automatically, yet, as the
% solutions of the linear systems might be computationally expensive, you
% must call the method 'expand' explicitly in order to update the
% expansions. Please note: If you call rp.expand in the command window,
% a list of all available quantities is displayed and you can select the
% one you are interested in with the mouse. The same holds for the plot
% function. 
rp.computeExpansion('Transmission')

% The plot of the error is updated automatically.

%% Adjust and refine the contour

rpc = RieszProjection(@f, parameters);

% As the poles are positioned on a tilted line, we try to make the contour 
% pass in the middle between two neighboring poles.
rpc.defineBg = [1.5-1.5i 6-1.5i];
rpc.minD = 0.0;
rpc.nPointsB = 128;
rpc.getContours;
fc1 = rpc.plot('ComplexPlane','fignumber',11);
rpc.computeExpansion('Transmission');
fc2 = rpc.plot('Error','fignumber',22);

rpc.computeExpansion('Transmission');

% Here, the rate of convergence is indeed much faster, and using the same
% number of integration points results in a lower error. This shows, that
% it is sometimes better to shift the contour instead of just doubling the
% number of integration points. Even if this requires a larger contour as
% in the given case.

%% Use eigenvalues from external source

% This time we assume that the eigenvalues are given and want to
% investigate the expansion over a larger energy interval.
% First, the quasi exact reference solution is obtained from the linearized
% problem. We request 30 eigenvalues centered at 6. 
ref = f([30,6],'Eigensystem');
evs = diag(ref.eigenvalues);

% We reset the results
rp.resetResults; rpc.resetResults;

% Now we can equip our instances with the quasi exact eigenvalues. The
% outer contour is kept and all eigenvalues inside the contour are
% selected automatically.
rp.poles = evs;
rpc.poles = evs;

% Now, we are interested in the frequency range from 2 to 13.
rpc.expansionPoints = 2:0.1:15;
rpc.referencePoints = rpc.expansionPoints(1:4:end);
rpc.defineBg = [1.5-1i 15-1i];

% We want to use the same poles and the same expansion points for the
% polygon and we must add the small contours again.
rp.expansionPoints = 2:0.1:15;
rp.referencePoints = rp.expansionPoints(1:4:end);
rp.defineBg = [1.5-1i 15-1i];

% We will start with a small number of integration points and refine the
% contours adaptively setting the property 'precision'. This property can
% be given as a scalar as below or as a vector with two entries. In the
% latter case, the first entry is used as an upper bound of the relative
% error and the second as an upper bound of the absolute value. For
% adaptive refinement, the minimum of absolute and relative error is used.
rpc.nPointsB = 60; rpc.nPoints = 8;
rp.nPointsB = 60; rp.nPoints = 8;
rp.precision = 1e-3;
rpc.precision = 1e-3;

% In this use case, we will compare the trapezoidal rule to the
% Gauss-Kronrod quadrature rule on a polygon.
rp.getContours('p',true) 
rp.minD = [0.6 2];
rpc.getContours('c');
rpc.minD = 0.6;

% Here, the second argument refers to the use of the Gauss-Kronrod
% quadrature rule.

%% 

%% Expand

rp.computeExpansion('Transmission');
rpc.computeExpansion('Transmission');

% If we have not closed the error plots, they will have been updated
% (otherwise we can just replot them) and we see that in this particular
% example, both quadrature rules yield results below the target precision
% after solving 120 linear systems along the background contour (for the
% (higher order quadrature rule 30 solutions cannot be reused after the
% refinement).
% Nevertheless, the exponential convergence of the trapezoidal rule, given
% circular or elliptic contours, will in practice often beat higher order
% quadrature rules.

%%
% The minimum of relative and absolute error, which has been used for the
% adaptive refinement can be plotted as well.
rp.plot('Error','errortype','m')

%% Function used to evaluate the scattering problem
function out = f(varargin)
% For custom solvers, you need a similar function that, depending on the
% input arguments, either solves a linear system (nargin==1) or applies a
% post process. In this particular case, we add a third option that returns
% the eigenspectrum of the linearized system. 
% The function is expected to return fields {{f11 f12 ...} {f21 ...} ...}
% if the input have been frequencies {[w11 ...] [w21 ...] ...} and if the
% input are fields {{f11 f21 ...} {f21 ...} ...} it is expected to return
% the corresponding target quantities {[q11 q12 ...] [q21 ...] ...}. If the
% quantity is quadratic in the fields the input is of the form:
% {{f111 f121 ...} {f112 f122 ...};{f211 f221 ...} {f212 f222 ...};...}. So
% the first index indicates the contour, the second the integration point
% on this contour, and the third distinguishes the normal field from the
% circfield. If the first input contains fields a second argument is
% expected that specifies the post process. 

persistent Ah Bh Ch jh s a;
if isempty(Ah)
    % Define sparse matrices used as a discrete representation of the problem.
    n = 300; L = pi*sqrt(2)/2; h=2*L/(n+1); V0 = 10;
    ndx1 = [2:n+2 1:n+1 1:n+2]; ndx2 = [1:n+1 2:n+2 1:n+2]; a = 1;
    Ah = sparse(ndx1,ndx2,h/6*[ones(1,2*n+2) 2 4*ones(1,n) 2]);
    Bh = sparse([1 n+2],[1 n+2],[1 1]);
    Ch = sparse(ndx1,ndx2,1/h*[-1*ones(1,2*n+2) 1 2*ones(1,n) 1]) - V0*Ah;
    % define physical source as a function of the frequency
    jh = @(w) [2i*w*exp(-1i*w*pi/sqrt(2)); zeros(n+1,1)]; 
    s = @(w,x) a*exp(1i*w*x);
end
% Get quasi exact solution to the linearized problem (nc: [# center]).
if nargin==2 && strcmpi('eigensystem',varargin{2})
    out = struct; nc = varargin{1}; if length(nc)==1, nc = [nc 0]; end
    E = speye(size(Ah)); O = sparse(size(Ah,1),size(Ah,1));
    A = [-1i*Bh, Ch; E,  O];  B = [Ah,  O;  O, E];
    [out.eigenvectors,out.eigenvalues]=eigs(A,B, nc(1), nc(2));
elseif iscell(varargin{1})
    if isempty(varargin{1})||ischar(varargin{1}{1}), out = {}; return; end
    out = cell(1,size(varargin{1},1+(nargin==1)));
    for it = 1:length(out)
        if nargin==1
            out{it} = solveScatteringProblems(varargin{1}{it});
        else
            out{it} = postProcess(varargin{1}(it,:));
        end
    end
end

    function out = solveScatteringProblems(ws)
        % Solve the inhomegeneous equation
        out = cell(1,numel(ws));
        for it2=1:numel(ws)
            w = ws(it2);
            out{it2} = {(w^2*Ah+1i*w*Bh-Ch)\(a*jh(w)) w};
        end
    end


    function out = postProcess(fs)
        if contains(lower(varargin{end}),'rightboundary')
            out = cellfun(@(x)x{1}(end),fs{1});
        elseif contains(lower(varargin{end}),'transmission')
            % The probability that the particle has passed the quantum well
            f1 = cellfun(@(x)x{1}(end),fs{1});
            if length(fs)==2
                % Because of the relation f(-w) = conj(f(conj(w)) the 
                % problem is not solved at negative but at complex
                % conjugate frequencies and therefore f(-w) we obtained
                % taking the complex conjugate.
                f2 = conj(cellfun(@(x)x{1}(end),fs{2}));
            else
                % In the case of the reference solution i.e. if the 
                % frequencies are real.
                f2 = conj(f1);
            end
            out = f1.*f2;
        elseif contains(lower(varargin{end}),'reflection')
            % The probability that the particle is scattered back.
            f1 = cellfun(@(x)x{1}(1),fs{1});
            if length(fs)==2
                % As in the case of transmission
                f2 = conj(cellfun(@(x)x{1}(1),fs{2}));
            else
                f2 = conj(f1);
            end
            w = cellfun(@(x)x{2},fs{1});
            out = (f1-s(w,-pi/sqrt(2))).*(f2-s(-w,-pi/sqrt(2)));
        end
    end
end
