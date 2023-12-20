classdef RieszProjection < handle
    %RIESZPROJECTION organization of projects for contour integration
    %   On the one hand, this class is designed to allow for the easy 
    %   construction of well-suited contours in the complex plain and 
    %   the application of different quadrature rules with adaptive
    %   refinement. On the other hand with the class '<a href="matlab: 
    %   help Scattering">Scattering</a>' an
    %   interface to the FEM solver JCMsuite is provided, that gives access
    %   to all the functionalities of this contour integration method.
    %   The aim is to perform contour integration methods for modal
    %   expansions and compute poles, zeros and corresponding derivatives, 
    %   making use of symmetries and parallelization for high efficiency.
    %   If the number of relevant eigenvalues is small (i.e., <= 20) the
    %   evaluation of the eigenvalues and the expansion can be done with a
    %   single contour enclosing the relevant region of the complex plane.
    %   Alternatively, expansions of arbitrary quantities based on the 
    %   electric field at any distance from the scatterer can be performed
    %   using quasi normal modes (QNMs) and corresponding expansion 
    %   coefficients. Normalized eigenmodes can be obtained with the 
    %   <a href="matlab: web(['https' ...
    %   '://jeos.springeropen.com/articles/10.1186/s41476-019-0098-z'], ...
    %   '-browser')">solver of JCMsuite</a>. Please refer to its <a 
    %   href="matlab: web(['https://www.docs.jcmwave.com/' ...
    %   'JCMsuite/html/ParameterReference/9184954f6f210148ef1d27247' ...
    %   'b8d7d47.html'],'-browser')">parameter reference</a> for 
    %   details.
    %   The idea to use a contour integral based approach dates back to a
    %   publication in 2018 by Zschiedrich et al. (<a
    %   href="matlab: web(['https://journals.aps.org/pra/abstract/' ...
    %   '10.1103/PhysRevA.98.043806'],'-browser')">PhysRevA</a>). The
    %   authors were inspired by Riesz projections (RPs), but the current
    %   implementation is rather based on contour integral methods for
    %   meromorphic functions. Not the total field is integrated but
    %   physical quantities are analytically continued to the complex plane
    %   and can be interpolated using a rational approximation based on the
    %   complex poles, which are the eigenvalues of the physical system.
    %   For further details on the expansion of quadratic functionals using
    %   RPs, including such related to the far field of the scatterer,  
    %   please refer to the article from 2020 by Binkowski et al. (<a
    %   href="matlab: web(['https://journals.aps.org/prb/abstract' ...
    %   '/10.1103/PhysRevB.102.035432'],'-browser')">PhysRevB</a>). 
    %   The implemented eigensolver is based on an algorithm reviewed
    %   by Austin et al. for integration contours based on roots of unity
    %   (<a href="matlab: web(
    %   'https://doi.org/10.1137/130931035','-browser')"
    %   >SIAM</a>). It is used by default if no eigenvalues are given.
    %   Please refer to the properties <a href="matlab: 
    %   help('RieszProjection.poles')">poles</a> and <a href="matlab: 
    %   help('RieszProjection.nMax')"
    %   >nMax</a>.
    %   
    %RIESZPROJECTION Properties:
    %   poles - Poles of f (see below) if they are known 
    %   selectedPoles - Define contributions
    %   expansionPoints - Expansion points on the real axis
    %   referencePoints - Frequencies for reference solutions
    %   nPoints - Number of points for all contours except the background
    %   nPointsB - Number of points for the background contour
    %   radius - Radius of circles enclosing the selected poles
    %   minD - Minimal distance of the background to the defining points
    %   contours - Discretized contours in the complex plane
    %   precision - Precision for adaptive refinement
    %   useEigenvectors - Expansion based on the normalized eigenmodes
    %   includeConjugatedPoles - Include the complex conjugate poles
    %   nMax - Upper abundand estimate of the number of eigenvalues
    %   cutOff - Select physcical poles based on the cut off frequencies
    %   nInterpolationNodes - Number of nodes for background interpolation
    %   defineBg - Complex numbers used to define the outer contour
    %   autoselect - Automatically select the poles inside the contour
    %   quantities - Quantities available for expansion
    %   f - Meromorphic function to be integrated along the contours
    %   fields - Cell array as returned by f
    %   derivedQuantities - Values at integration points, e.g., energies
    %   reference - Struct containing reference solutions
    %   contributions - Struct holding expansion data, e.g., residues
    %   zeros - Stuct holding zeros
    %   derivatives - Partial derivatives if available
    %   nModes - The number of considered modes
    %
    %RIESZPROJECTION Methods:
    %   RIESZPROJECTION - Construct an instance with given properties
    %   selectPoles - Select poles interactively
    %   getContours - Construct the contours for integration
    %   getWeights - Weights for quadrature
    %   computeExpansion - Expand a selected quantity
    %   computeReference - Evaluate reference solutions at real omega
    %   computeError - Evaluate errors of the integrals and the expansion
    %   computePoles - Compute poles, i.e., eigenvalues of the system
    %   computeZeros - Compute zeros of the specified quantity
    %   plot - Plot results and errors
    %
    %RIESZPROJECTION static Methods:
    %   minimalCircle - Compute minimal area enclosing circle
    %   minimalEllipse - Compute minimal area enclosing ellipse
    %
    %DEPENDENCIES:
    %   The class 'RieszProjection' can be used with any solver that can
    %   provide the analytic continuation of the observable of interest to
    %   the complex plane. An interface in the form of a callable must be
    %   passed to the constructor of this class. Furthermore, information
    %   about the quantities that are available for expansions must be
    %   provided by setting the property <a href="matlab: 
    %   help('RieszProjection.quantities')">quantities</a>. 
    %   There are a few lines, which are specific for the interface
    %   <a href="matlab: help('Scattering')"
    %   >Scattering</a> to the FEM solver JCMsuite. For use with this
    %   interface, third-party support has to be set up.
    %   Please refer to the documentation of the <a href="matlab:
    %   web(['https://docs.jcmwave.'com/JCMsuite/html/MatlabInterface' ...
    %   '/index.html'],'-browser')">matlab interface</a> to set up
    %   the third-party support and to the homepage of <a href="matlab:
    %   web('https://jcmwave.com/jcmsuite/jcmsuite-documentation', ...
    %   '-browser')">JCMwave</a> for
    %   installation instructions.
    %
    %   see also Scattering

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: April-2023
      
    properties (Dependent)
        % poles - Poles of the system if known. If no poles are given,
        % the eigenvalue problem is solved using a contour integral
        % algorithm. This requires an estimated number of eigenvalues
        % expected to be found inside the contour. This can be set with the
        % property <a href="matlab: help( ...
        % 'RieszProjection/nMax')">nMax</a> and should rather be chosen too
        % large than too small. In the case of underestimating the number 
        % of eigenvalues, the solutions inside the contour are false.
        poles (1,:) double;
        
        % selectedPoles - Grouped poles define individual contributions.
        % This cell array contains vectors of integers which are indices to
        % the property 'poles'. If you do not want to select any, set this 
        % property to {0}. Each entry, i.e, selectedPoles{k} can either be
        % a scalar or a vector. If it is a vector, the contributions of the
        % corresponding poles are summed up.
        selectedPoles (1,:) cell;
        
        % expansionsPoints - Calling the method 'computeExpansion' will
        % return the expansion of the quantity of interest at the expansion
        % points defined here, unless otherwise specified. 
        expansionPoints (:,1) double;
        
        % referencePoints - If not empty the scattering problems at these
        % points are solved and used to estimate the error of the
        % expansion.
        referencePoints (:,1) double;
        
        % nPoints - Number of points on the contours around the selected 
        % poles. It can be a scalar or an array of the same length as the
        % property 'selectedPoles'. The default is zero, i.e., usually
        % these contours are not required, and solving the eigenvalue
        % problem and evaluating the expansion can be based exclusively on
        % a single contour.
        nPoints (1,:) double {mustBeNonnegative,mustBeInteger};
       
        % The number of points defining the background contour. If this
        % property is set to zero, e.g., if normalized eigenmodes are
        % available from an external source, the difference between the
        % background and some reference points is interpolated using
        % splines.
        nPointsB (1,1) double ...
            {mustBeNonnegative, mustBeInteger};
        
        % radius - Radius of the circles enclosing the selected poles. If
        % more than one pole must be enclosed it defines the minimal
        % distance of the contour to any of the poles
        radius (1,:) double {mustBeNonnegative};
        
        % minD - Approximate minimal distance of the background contour to 
        % enclosed poles. If the contour is an ellipse it can be a
        % two-element vector. The first value is then added to the
        % semi-major axis and the second to the semi-minor axis. In the
        % case of a polygon the first entry results in a shift in x-
        % direction and the second entry in a shift in y-direction.
        minD (1,:) double;
        
        % contours - Contours associated with selected poles. The
        % background contour RP.contours{end} is the contour enclosing all
        % poles of interest and all expansion points if given.
        contours (1,:) cell;
        
        % precision - This property can be a scalar or a vector with two
        % elements. In the latter case, the first value is used as the upper
        % bound of the absolute error and the second one as upper bound of
        % the relative error. In the former case, both bounds have the same
        % value. This property is used for an adaptive refinement of the
        % integration contours. For absolute and relative error control
        % you can set the first or the second value to zero respectively.
        % If the precision is Inf (default), no refinements are made. 
        precision (1,2) double;
        
        % useEigenvectors - If this property is set to true the normalized
        % eigenvectors are used together with the analytic expansion
        % coefficients if given.
        useEigenvectors (1,1) logical;
        
        % includeConjugatedPoles - If set to true, the poles of f(-w) are
        % taken into account which is necessary for the expansion of
        % quadratic forms such as the electric field energy.
        includeConjugatedPoles (1,1) logical;
        
        % nMax - If the eigenvalues are not given, a guess
        % of the expected number of eigenvalues is needed. It should better
        % be too large, than too small. If 3 eigenvalues are contained in
        % the contour nMax = 6 might be a good choice, if you expect 10
        % eigenvalues, nMax = 15 might be adequate. Keep in mind,
        % that for quadratic quantities all poles in the lower part of the
        % complex plane will have a counterpart above the real line. If you
        % call the method 'computePoles' with a quadratic
        % quantity the dimension of the Hankel matrix will be 2*nMax.
        % Eigenvalues whose normalized residues are smaller than a 
        % specified cutoff are discarded. If this property is set, the 
        % property 'poles' is reset to an empty list. 
        nMax (1,1) double {mustBePositive,mustBeInteger};
        
        % cutOff - Due to noise, usually nMax eigenvalues based on higher
        % moments are evaluated, but not all of them correspond to
        % eigenvalues of the physical system. The residues can be used to
        % measure the coupling of the corresponding quantity to the given
        % source. The residues are normalized with respect to the largest
        % one and all eigenvalues with a normalized residue larger than the
        % cut off are discarded. 
        cutOff (1,1) double {mustBePositive};
        
        % nInterpolationNodes - If nPointsB = 0, the background contribution 
        % is estimated by interpolation using cubic splines.
        nInterpolationNodes (1,1) double {mustBeNonnegative,mustBeInteger};
        
        % defineBg - If not empty, these points are used instead of the
        % selected poles in order to define the background contour.
        defineBg (:,:) double;
        
        % autoselect - If true all eigenvalues contained in the contour are
        % selected automatically. The default is true if the background
        % contour is used for the residual contribution and false
        % otherwise.
        autoselect (1,1) logical; 
        
        % quantities - Quantities that are available for expansion.
        % Working with JCMsuite and the class Scattering, this property is
        % set automatically. Otherwise, its fields are expected to be
        % structs with additional information. The field 'quadratic' is
        % mandatory for each quantity. Furthermore the optional field
        % 'evaluate' may contain a function handle whose first argument is
        % the expansion frequency and the following arguments are defined
        % with the field 'parents'. If the latter is missing the parent is
        % assumed to be the quantity itself. E.g.:
        % p = {'Radiation','DipoleEmission'};
        % f = @(~,rad,dpe)(rad./sum(dpe,3);
        % quantities.eta = struct('parent',p,'evaluate',f,'quadratic',true)
        quantities (1,1) struct {mustBeValidQuantities};
    end
    
    properties (SetAccess=private)
        % f - Function of complex frequencies w to be integrated along the
        % contours. You can provide f either as a function handle or a
        % callable class which takes a cell array of complex doubles with 
        % shape=(1,:).
        % Working with JCMsuite you should use an instance of the class
        % Scattering which returns a cell array containing cells of
        % character arrays. This property is set by the constructor.
        % f({}) must return {} to pass the simple test mustTakeCell.
        f {mustTakeCell};
        
        % fields - Cell array as returned by f
        fields (1,:) cell;
        
        % derivedQuantities - Quantities derived from the fields, e.g., 
        % energies, corresponding to function values of the meromorphic
        % function at the integration points.
        derivedQuantities (1,1) struct;
        
        % reference - Struct containing reference solutions in order to get
        % the reference solutions at the frequencies w0 you must call the
        % method <a href="matlab:
        % help('RieszProjection/computeReference')"
        % >computeReference</a>. 
        reference (1,1) struct;
        
        % eigenmodes - Data required for the modal contributios is stored
        % in this struct. It has fields corresponding to the expanded
        % quantities.
        contributions (1,1) struct;
        
        % zeros - Struct containing the zeros which are evaluated calling
        % the method 'computeZeros'. The fieldnames correspond to quantities.
        zeros (1,1) struct; 
        
        % derivatives - Derivatives of poles and zeros with respect to 
        % design parameters. They are evaluated if the linear quantity is 
        % given with derivatives.
        derivatives (1,1) struct;
        
        % nModes - The number of modes contributing to the modal picture
        % This property is updated at each call of evaluate.
        nModes = 0;
    end
    
    properties (Constant, Hidden)
        eps0 = 8.85418781762039e-12; % permittivity in vacuum
        mu0 = 4e-7*pi; % permeability in vacuum
        c0 = 299792458; % speed of light
        gn = [7 15 25];
        % list of properties
        % nodes and weights for common Gauss-Kronrod quadrature rules
        % https://www.advanpix.com/2011/11/07/
        % gauss-kronrod-quadrature-nodes-weights/
        nodes15 = [ ...
            0.2077849550078985; 0.4058451513773972; ...
            0.5860872354676911; 0.7415311855993944; ...
            0.8648644233597691; 0.9491079123427585; ...
            0.9914553711208126];
        weights15 = [ ...
            0.2094821410847278, 0.2044329400752989, ...
            0.1903505780647854, 0.1690047266392679, ...
            0.1406532597155259, 0.1047900103222502, ...
            0.06309209262997855, 0.02293532201052922];
        weights7g = [ ...
            0.4179591836734694, 0, 0.3818300505051189, 0, ...
            0.2797053914892767, 0, 0.1294849661688697, 0];
        nodes31 = [ ...
            0.1011420669187175; 0.2011940939974345; ...
            0.2991800071531688; 0.3941513470775634; ...
            0.4850818636402397; 0.5709721726085388; ...
            0.6509967412974170; 0.7244177313601700; ...
            0.7904185014424659; 0.8482065834104272; ...
            0.8972645323440819; 0.9372733924007059; ...
            0.9677390756791391; 0.9879925180204854; ...
            0.9980022986933971];
        weights31 = [ ...
            0.1013300070147915, 0.1076984552387560, ...
            0.09917359872179196, 0.09664272698362368, ...
            0.09312659817082532, 0.08856444305621177, ...
            0.08308050282313302, 0.07684968075772038, ...
            0.06985412131872826, 0.06200956780067064, ...
            0.05348152469092809, 0.04458975132476488, ...
            0.03534636079137585, 0.02546084732671532, ...
            0.01500794732931612, 0.005377479872923349];
        weights15g = [ ...
            0.2025782419255613, 0, 0.1984314853271116, 0, ...
            0.1861610000155622, 0, 0.1662692058169939, 0, ...
            0.1395706779261543, 0, 0.1071592204671719, 0, ...
            0.07036604748810812, 0, 0.03075324199611727, 0];
        nodes51 = [ ...
            0.06154448300568508; 0.1228646926107104; ...
            0.1837189394210489; 0.2438668837209884; ...
            0.3030895389311078; 0.3611723058093878; ...
            0.4178853821930377; 0.4730027314457150; ...
            0.5263252843347192; 0.5776629302412230; ...
            0.6268100990103174; 0.6735663684734684; ...
            0.7177664068130844; 0.7592592630373576; ...
            0.7978737979985001; 0.8334426287608340; ...
            0.8658470652932756; 0.8949919978782754; ...
            0.9207471152817016; 0.9429745712289743; ...
            0.9616149864258425; 0.9766639214595175; ...
            0.9880357945340772; 0.9955569697904981; ...
            0.9992621049926098];
        weights51 = [ ...
            0.06158081806783294, 0.06147118987142532, ...
            0.06112850971705305, 0.06053945537604586, ...
            0.05972034032417406, 0.05868968002239421, ...
            0.05743711636156783, 0.05595081122041232, ...
            0.05425112988854549, 0.05236288580640748, ...
            0.05027767908071567, 0.04798253713883671, ...
            0.04550291304992179, 0.04287284502017005, ...
            0.04008382550403238, 0.03711627148341554, ...
            0.03400213027432934, 0.03079230016738749, ...
            0.02747531758785174, 0.02400994560695321, ...
            0.02043537114588284, 0.01684781770912830, ...
            0.01323622919557167, 0.009473973386174152, ...
            0.005561932135356714, 0.001987383892330316];
        weights25g = [ ...
            0.1231760537267155, 0, 0.1222424429903100, 0, ...
            0.1194557635357848, 0, 0.1148582591457116, 0, ...
            0.1085196244742637, 0, 0.1005359490670506, 0, ...
            0.09102826198296365, 0, 0.08014070033500102, 0, ...
            0.06803833381235692, 0, 0.05490469597583519, 0, ...
            0.04093915670130631, 0, 0.02635498661503214, 0, ...
            0.01139379850102629, 0];
    end
    
    % Get and set methods for the dependent properties.
    % If a dependent property is set to a new value, contours will be
    % reevaluated and function values which need to be
    % recalculated are marked in the cell array undefinedValues.
    methods
        function set.nPoints(RP, nPoints)
            if isequal(nPoints, RP.nPoints), return; end
            if ~isempty(RP.selectedPoles)
                n = length(RP.selectedPoles);
                if isscalar(RP.radius_)
                    RP.radius_ = repmat(RP.radius_,1,n);
                end
                if isscalar(nPoints), nPoints = repmat(nPoints,1,n);
                elseif length(nPoints)~=length(RP.selectedPoles)
                    error(['To set the property ''nPoints'' provide a '...
                        'scalar or an array which matches the length '...
                        'of the property ''selectedPoles'' (%d)'], ...
                       	length(RP.selectedPoles));
                end
            end
            RP.nPoints_ = nPoints; 
            RP.getContours(true);
        end
        function value = get.nPoints(RP)
            if all(RP.nPoints_==RP.nPoints_(1))
                value = RP.nPoints_(1);
            else
                value = RP.nPoints_;
            end
        end
        
        function set.radius(RP, radius)
            if isequal(radius,RP.radius), return; end
            if ~isempty(RP.selectedPoles)
                if isscalar(radius)
                    radius = repmat(radius,1,length(RP.selectedPoles));
                elseif length(radius)~=length(RP.selectedPoles)
                    error(['To set the property ''radius'' provide a '...
                        'scalar or an array which matches the length '...
                        'of the property ''selectedPoles'' (%d)'], ...
                        length(RP.selectedPoles));
                end
            end
            RP.radius_ = radius; 
            RP.getContours(true);
        end
        function value = get.radius(RP)
            if all(RP.radius_==RP.radius_(1))
                value = RP.radius_(1);
            else
                value = RP.radius_;
            end
        end
        
        function set.minD(RP, minD)
            if all(minD==RP.minD), return; end
            if isscalar(minD), minD = [minD minD]; end
            RP.minD_ = minD;
            RP.getContours(false);
        end
        function value = get.minD(RP)
            value = RP.minD_;
        end
        
        function set.nPointsB(RP, nPointsB)
            if ~isempty(RP.contours)
                nP = numel(RP.contours{end});
            else
                nP = RP.nPointsB_;
            end
            if nPointsB==nP, return; end
            RP.nPointsB_ = nPointsB;
            if nPointsB, RP.nInterpolationNodes = 0; end
            % only the number of nodes will change but not the shape
            RP.getContours(true); 
        end
        function value = get.nPointsB(RP)
            value = RP.nPointsB_;
        end
        
        function set.expansionPoints(RP, expansionPoints)
            if isequal(expansionPoints,RP.expansionPoints_), return; end
            expansionPoints = sort(expansionPoints);
            try ps1 = RP.expansionPoints_([1 end]); catch, ps1=[]; end
            try ps2 = expansionPoints([1 end]); catch, ps2 = []; end
            RP.expansionPoints_ = expansionPoints(:);
            tolerance = -RP.minD + min(1e-2,RP.precision(1))*RP.minD;
            noBg = ~isempty(RP.fields) && all(inside(RP,ps2,tolerance));
            if noBg || isempty(RP.contours) || isequal(ps1,ps2)
                updatePlots(RP)
            else
                RP.getContours(false);
            end
        end
        function value = get.expansionPoints(RP)
            value = RP.expansionPoints_;
        end
        
        function  set.referencePoints(RP, referencePoints)
            if isequal(referencePoints,RP.referencePoints_), return; end
            RP.referencePoints_ = sort(referencePoints);
            RP.up2date(5) = false;
        end
        function value = get.referencePoints(RP)
            value = RP.referencePoints_;
        end
        
        function set.poles(RP, poles)
            poles = sort(poles,'ComparisonMethod','real');
            if isequal(RP.poles, poles), return; end
            RP.poles_ = poles; 
            RP.up2date(2:3) = ~isempty(poles);
            RP.contributions = struct;
            if isempty(RP.contours)
                updatePlots(RP)
                return;
            end
            sP = find(inside(RP,poles));
            RP.selectedPoles = num2cell(sP);
        end
        function value = get.poles(RP)
            value = RP.poles_;
        end
        
        function set.selectedPoles(RP, selectedPoles)
            if isequal(selectedPoles, RP.selectedPoles), return; end
            if ~isempty(selectedPoles)
                selectedPoles = cellfun(@(x){x(:).'},selectedPoles);
            end
            resetSelectedPoles(RP,selectedPoles);
            if ~isempty(RP.selectedPoles)
                ndx = [selectedPoles{:}];
                tolerance = min(1e-2,RP.precision(1))*RP.minD;
                inC = inside(RP,RP.poles(ndx),tolerance);
                if all(inC) && sum(inC)==sum(inside(RP,RP.poles))
                    RP.getContours(true);
                else
                    RP.getContours(false);
                end
            else
                RP.up2date(2) = true;
            end
        end
        function value = get.selectedPoles(RP)
            value = RP.selectedPoles_;
        end
        
        function set.contours(RP, contours)
            if isempty(RP.contours_) || isequal(contours, RP.contours_)
                RP.contours_ = contours;
                if isempty(RP.selectedPoles) && ~isempty(RP.poles)
                    RP.selectedPoles = num2cell(find(inside(RP,RP.poles)));
                end
                updatePlots(RP)
                return;
            end
            RP.contours_ = contours;
            if RP.autoselect
                n = length(RP.selectedPoles);
                resetSelectedPoles(RP,{}); 
                if RP.nPoints && n~=length(RP.selectedPoles), return; end
            end
            updatePlots(RP)
            if isempty(RP.contours_f), return; end
            RP.undefinedValues = cell(1,length(contours));
            RP.definedValues = cell(1,length(contours));
            c = cellfun(@(x){x(:)},RP.contours_f); c = cat(1,c{:});
            c = cat(2,real(c),imag(c));
            for it = 1:length(contours)
                c_it = contours{it}(:);
                if isempty(c_it), continue; end
                c_it = cat(2,real(c_it),imag(c_it));
                [def,loc] = ismembertol(c_it, c, 'ByRows', true);
                RP.undefinedValues{it}=~def.';
                RP.definedValues{it} = loc(def).';
            end
            RP.up2date(1) = false;
            % only update the eigenvalues if they are not given from an
            % external source and if the precision is expected to increase
            decreased = sum(cellfun(@numel,contours))<size(c,1);
            RP.up2date(2) = RP.up2date(2) && decreased && ~RP.nPoints;
        end
        function value = get.contours(RP)
            value = RP.contours_;
        end
        
        function set.precision(RP, p)
            this = 'RieszProjection.precision';
            if isscalar(p)
                if ~isfloat(p) && isreal(p) && (p > 0)
                    error('Invalid value for %s', this);
                end
                RP.precision_ = [p p];
            elseif numel(p)==2
                pabs = p(1); prel = p(2);
                if ~isfloat(pabs) && isreal(pabs) && (pabs > 0)
                    error('Invalid absolute tolerance for %s', this);
                elseif ~isfloat(prel) && isreal(prel) && (prel > 0)
                    error('Invalid relative tolerance for %s', this);
                end
                RP.precision_ = [pabs prel];
            end
        end
        function value = get.precision(RP)
            if RP.precision_(1)==RP.precision_(2)
                value = RP.precision_(1);
            else
                value = RP.precision_;
            end
        end
        
        function set.includeConjugatedPoles(RP,iCP)
            if iCP==RP.includeConjugatedPoles, return; end
            RP.includeConjugatedPoles_ = iCP;
            RP.getContours(false);
            aQs = fieldnames(RP.quantities);
            if ~RP.includeConjugatedPoles_
                aQs = aQs(structfun(@(x)~x.quadratic,RP.quantities));
            end
            RP.availableQuantities = aQs;
        end
        function value = get.includeConjugatedPoles(RP)
            value = RP.includeConjugatedPoles_;
        end
        
        function set.useEigenvectors(RP,uEV)
            if uEV==RP.useEigenvectors, return; end
            RP.useEigenvectors_ = uEV;
            if ~uEV && RP.nPoints_(1)
                uV = cellfun(@(x){true(size(x))},RP.contours(1:end-1));
                if isempty(RP.undefinedValues)
                    sz = [1 max(length(RP.contours),1)];
                    RP.undefinedValues = cell(sz);
                    RP.definedValues = cell(sz);
                end
                RP.undefinedValues = [uV RP.undefinedValues(end)];
                RP.definedValues(1:end-1) = cell(size(uV));
                RP.up2date(1) = false;
            end
        end
        function value = get.useEigenvectors(RP)
            value = RP.useEigenvectors_;
        end
        
        function set.nMax(RP,nEv)
            if nEv==RP.nMax_, return; end
            RP.nMax_ = nEv;
            if ~isempty(RP.poles)
                RP.up2date(2) = false;
            end
        end
        function value = get.nMax(RP)
            value = RP.nMax_; 
        end
        
        function set.cutOff(RP,cO)
            if cO==RP.cutOff_, return; end
            RP.cutOff_ = cO;
            if ~isempty(RP.poles)
                RP.up2date(2) = false;
            end
        end
        function value = get.cutOff(RP)
            value = RP.cutOff_; 
        end
        
        function set.nInterpolationNodes(RP,nIN)
            if nIN==RP.nInterpolationNodes_, return; end
            RP.nInterpolationNodes_ = nIN;
            c = RP.contributions; ns = fieldnames(c).';
            for n = ns, c.(n{1}) = rmfield(c.(n{1}),'interpolation'); end
            RP.contributions = c;
            if nIN
                RP.nPointsB = 0;
                RP.autoselect = false;
            end
            RP.up2date(5) = false;
        end
        function value = get.nInterpolationNodes(RP)
            value = RP.nInterpolationNodes_;
        end
        
        function set.defineBg(RP,defineBg)
            if isequal(defineBg(:),RP.defineBg), return; end
            RP.defineBg_ = defineBg(:);
            RP.getContours(false);
        end
        function value = get.defineBg(RP)
            value = RP.defineBg_;
        end
        
        function set.autoselect(RP,autoselect)
            if isequal(autoselect,RP.autoselect_), return; end
            RP.autoselect_ = autoselect;
            resetSelectedPoles(RP,{});
        end
        function value = get.autoselect(RP)
            value = RP.autoselect_;
        end
        
        function set.quantities(RP,qs)
            ch = struct;
            for n = fieldnames(qs).'
                if ~isfield(qs.(n{1}),'hidden')
                    qs.(n{1}).hidden = false;
                end
                if ~isfield(qs.(n{1}),'m')
                    qs.(n{1}).m = 1; % q is assumed to be scalar
                end
                if isfield(qs.(n{1}),'parents')
                    ch.(qs.(n{1}).parents{1}).(n{1}) = [];
                end
            end
            for n = fieldnames(ch).'
                qs.(n{1}).children = fieldnames(ch.(n{1})).';
            end
            RP.quantities_ = qs;
            aQs = fieldnames(qs);
            if ~RP.includeConjugatedPoles
                aQs = aQs(structfun(@(x)~x.quadratic,qs));
            end
            RP.availableQuantities = aQs;
        end
        function value = get.quantities(RP)
            value = RP.quantities_;
        end
    end
    
    properties (Access=private)
        % See dependent properties
        nPoints_ = 0;
        nPointsB_ = 64;
        radius_ = 1e13;
        minD_ = 2e13;
        expansionPoints_;
        referencePoints_;
        poles_;
        selectedPoles_ = {};
        contours_ = {};
        precision_ = Inf(1,2);
        useEigenvectors_ = false;
        includeConjugatedPoles_ = false;
        nMax_ = 12;
        cutOff_ = 1e-2;
        nInterpolationNodes_ = 0;
        defineBg_;
        autoselect_ = true; 
        quantities_ = struct;
        name_; % the name of the instance for hrefs
        
        % names of the quantities available for expansion
        availableQuantities (1,:) cell;
        
        % contours_f - The contours corresponding to the current fields
        % which are updated at each call of evaluate.
        contours_f;
        
        contourSettings (1,1) struct = struct;
        plotSettings (1,1) struct = struct('fignumber',struct);
        % 1 - integration points
        % 2 - eigensystem (only update if accuracy is increased)
        % 3 - eigenvalues from external source
        % 4 - eigenvectors from external source
        % 5 - plots
        up2date (1,5) logical = false(1,5);

        % Values which are not defined in RP.fields are marked and
        % locations of already existing values are saved in definedValues
        undefinedValues (1,:) cell;
        definedValues (1,:) cell;
        
        % shape - saves information about the shape of the background
        % shape{1} - circle: [c r], ellipse: [0 1], polygon: vertices
        % shape{2} - 'c', 'e' or 'p'
        % shape{3} - order of Gauss-Kronrod quadrature rule
        % shape{4} - ellipse parameters or vertices without splits (pol.)
        shape (1,:) cell;
        
        % weights - Defined by the contours and needed for the quadrature
        % rule. For a contour gamma those are the derivative of gamma at
        % the corresponding point times some value defined by the selected
        % quadrature rule, e.g., the spacing h between the integration
        % points in case of the trapezoidal rule.
        weights (1,:) cell;
    end
    
    methods
        function RP=RieszProjection(f, varargin)
            % Construct an instance of the class
            %    rp = RIESZPROJECTION(f) Constructs an instance with
            %       default properties based on evaluations of f at
            %       complex frequencies.
            %    rp = RIESZPROJECTION(f,parameters) Constructs an instance
            %       whose properties are set according to the fields of
            %       the struct parameters. 
            %    rp = RIESZPROJECTION(f,PARAM1,VAL1,PARAM2,VAL2,...)
            %       Instead of a struct you can set the properties of this
            %       class using name value pairs where the names must be
            %       unique case insensitive partial matches of the property
            %       names
            %        
            % The function f defines the computationally expensive part of 
            % the function to be integrated. The analytic part is given by
            % the coefficients a.
            
            RP.f = f;
            RP.reference = struct('frequencies',zeros(0,2),'fields',[]);
            if isa(f,'Scattering')
                if ~isempty(RP.f.evsDir)
                    RP.up2date(3:4) = true;
                    RP.poles = unique(RP.f.eigenvalues);
                end
                f.quantitiesForExpansion(RP); % sets the quantities
            end
            if nargin==1, return; end
            ps = properties(RP);
            ps_ = ps(1:find(strcmp(ps,'f'))-1);
            if isstruct(varargin{1})
                fs = fieldnames(varargin{1});
                vs = struct2cell(varargin{1});
            elseif length(varargin)>1
                fs = varargin{1:2:end};
                vs = varargin{2:2:end};
            end
            for it = 1:length(fs)
                p = validatestring(fs{it},ps_);
                RP.(p) = vs{it};
            end
        end
        
        % selectPoles - Select poles interactively
        selectPoles(RP, radius)
        
        % getContours - Get the contours according to the selected poles    
        getContours(RP, varargin);
        
        % getWeights - Get the quadrature weights
        weights = getWeights(RP,c,p);
        
        % plot - calls plotting functions for contours and quantities
        varargout = plot(RP, varargin)
            
        % computeExpansion - expand a given quantity
        varargout = computeExpansion(RP, varargin)
            
        % Computes the reference solution on the real axis
        varargout = computeReference(RP, quantity, w0, keys)
        
        % error - the maximum error along all expansion points is returned
        [err,x] = error(RP, quantity, index, kind)
        
        % computePoles - compute the poles (eigenvalues) of the system
        varargout = computePoles(RP,varargin)
        
        % computeZeros - compute the zeros of a specified quantity
        varargout = computeZeros(RP,varargin)
        
        function resetResults(RP)
            % RESETRESULTS reset all results
            RP.selectedPoles_ = {};
            RP.contours = {};
            RP.undefinedValues = {};
            RP.definedValues = {};
            RP.up2date(1) = false;
            RP.up2date(2) = RP.up2date(3);
            RP.contributions = struct;
            RP.derivatives = struct;
            RP.zeros = struct;
        end
        
        % getDir - directory of integration point or reference solution by 
        % contour index and the frequency of the point. The first argument 
        % is the index. If there are n modal contributions the index n+1 
        % refers to the background contour and the index n+2 to the
        % reference solution. The path to the scattering problem with the
        % frequency closest to the given one is returned.
        % This method is JCMsuite specific.
        varargout = getDir(RP,index,omega)
    end
    
    methods (Static)
        % Get the circle with minimal area enclosing all points.
        [c, r] = minimalCircle(points);
        % Get the ellipse with minimal area enclosing all points.
        E = minimalEllipse(points);
    end
    
    methods (Access=private)
        % Evaluate all function values on the contours.
        evaluate(RP)
        
        % Evaluate derived quantities. This method should only be called by
        % the method 'evaluate'
        updateDerivedQuantities(RP, d)
        
        % Evaluate modal contributions
        out = computeContributions(RP, varargin)
        
        % compute the eigenvalues
        varargout = locatePoles(RP,index,cutoff,nmax,q)
        
        % apply quadrature rule
        [res, convData] = quad(RP, values, vars, conv, it)
        
        % use contour integral to compute zeros and poles
        [ps,rs,ds,p] = zerosPoles(RP,poles,varargin)
    end
end

% Minimal test for f
function mustTakeCell(f)
if isempty(f), return; end
try a=f({}); catch, a = 1; end
if iscell(a) && ~isempty(a)
    error(['f should be callable taking cell arrays of doubles with '...
        'shape=(1,:) containing complex frequencies and return a'...
        ' cell array of the same size containing the results of the'...
        ' linear systems. At a second call with the previous output as'...
        ' input plus a character array defining a quantity, this'...
        ' quantity should be derived from the solutions of the' ...
        ' linear systems. f({}) must return {}']);
end
end

% test for quantities
function mustBeValidQuantities(quantities)
message = 'The validation of the quantity %s failed: %s';
qs = fieldnames(quantities);
for it = 1:length(qs)
    q = quantities.(qs{it});
    if ~isfield(q,'quadratic')
        error(message,qs{it},'Field ''quadratic'' missing')
    elseif ~islogical(q.quadratic)
        error(message,qs{it},'The field ''quadratic'' must be logical')
    end
    if isfield(q,'evaluate') && ~isa(q.evaluate,'function_handle')
        msg = 'The field ''evaluate'' must be a function_handle';
        error(message,qs{it},msg)
    end
    if isfield(q,'parents') 
        if ~all(ismember(q.parents,qs))
            msg = 'The ''parents'' must be valid quantities';
            error(message,qs{it},msg);
        elseif ismember(qs{it},q.parents)
            msg = 'A quantity must not be its own parent';
            error(message,qs{it},msg)
        elseif ~isfield(q,'evaluate')
            msg = 'The field ''parents'' requires the field ''evaluate''';
            error(message,qs{it},msg)
        end
    end
    if isfield(q,'children') 
        if ~all(ismember(q.children,qs))
            msg = 'The ''children'' must be valid quantities';
            error(message,qs{it},msg);
        elseif ismember(qs{it},q.children)
            msg = 'A quantity must not be its own child';
            error(message,qs{it},msg)
        end
    end
end
end
