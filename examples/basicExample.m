% Minimal example to get familiar with Riesz projection expansions.
% Please keep in mind that you can always use the help function of matlab
% to get detailed descriptions of the properties and methods you are
% interested in. E.g., execute help RieszProjection to get a list of its
% methods. You can follow the corresponding link to get more details about
% a particular method.

% clear the workspace
clearvars
clc
range('default');

%% Set up the Matlab interface to JCMsuite
% Activate the third party support and start a daemon. Please refer to the
% README for details.

% Provide the path to the installation directory of JCMsuite and add it to
% the Matlab path. The version should be 4.4.0 or higher.
jcm_root = '~/Software/JCMsuite.5.4.3';
addpath(fullfile(jcm_root, 'ThirdPartySupport', 'Matlab'));

% The directories '@RieszProjection', and '@Scattering' must be accessible. 
% Here, it is assumed that they are contained in the parent directory of
% this file.
addpath('../src');

% Start a daemon that handles the jobs submitted to jcmwave_solve.
options = struct('Hostname', 'localhost', ...
                'Multiplicity',8, ...
                'NThreads',2 ... 
                );     
            
% Shutdown a possibly running daemon and register a new computer ressource.
jcmwave_daemon_shutdown;
jcmwave_daemon_add_workstation(options);

% The working directory for RPExpand can be chosen independently from the
% location of the project files for JCMsuite. 
wDir = 'DiamondNanodisk';

% This directory will contain the solution of the eigenvalue problem.
resDir = [wDir filesep 'resonance' filesep 'project_results'];

%% Construct instance of Scattering
% This object takes care of solving the scattering problems for different
% frequencies while all other parameters remain constant. You can add all
% parameters you might need for parameter substitution to the structure
% 'keys' which is passed to the constructor together with the path to the
% project file.

% The struct 'keys' must have the field 'finiteElementDegree'. If you want
% to simulate a point source and want to expand its dipole emission Gamma,
% you must define position and strength in order to use the expression
% Gamma(w) = -real(E(w,position).'*strength) which assumes the strength to
% be real. Alternatively, you can expand  Gamma using the electromagnetic
% field energy flux, which is quadratic in the electric field strength E.

% Set the finite element degree for the scattering problem defined in
% 'example/project.jcmpt'.
keys.finiteElementDegree = 4;

% In this example, the dipole emission Gamma will be expanded. As stated
% above, in order to do this based on a point evaluation of the electric
% field, you have to define a strength vector, i.e., amplitude and
% polarization of the point source and its position within the 
% computational domain. Both parameters are used in the source file 
% 'example/sources.jcmt' and the latter additionally defines the position 
% where E is evaluated.
keys.position = [0 145e-9 0];
keys.strength = [1 0 0];

% Path to the project file. The line of code that follows this comment
% assumes, that it is stored in the directory 'scattering' located in the
% working directory 'wDir' defined above. 
projectFile = [wDir filesep 'scattering' filesep 'project.jcmpt'];

% Here, the constructor of the class 'Scattering' is called. All results are
% written to the directory 'wDir', which is passed to the constructor.
sc = Scattering(projectFile, keys, wDir);
sc.addQuantity('DipoleEmission')
x = (-350:5:350)*1e-9; y = (-250:5:450)*1e-9; z = 0;
export_keys = struct('GridPointsX',x,'GridPointsY',y,'GridPointsZ',z);
sc.addQuantity('CartesianExport',export_keys);

%% Construct an instance of RieszProjection
% You can pass a struct with options to the constructor which will be
% assigned to the corresponding properties. Alternatively, you can set them
% up after construction. For a complete list of public properties, type
% help RieszProjection into the command window.

% In this example, we are only interested in the expansion of the dipole
% emission given as an expression that is linear in the electric field.
% Symmetry with respect to the real axis is therefore not required and the
% property 'includeConjugatedPoles' can be set to false. If you want to 
% expand quantities that are quadratic in E, such as the electric field 
% energy, this property must be true. For further details refer to the 
% example quadraticForms.m.
parameters.includeConjugatedPoles = false;

% The values will be added to the axes of the ellipse or to the radius
% of the circle. In the latter case, it defines the minimal distance of the
% background contour to the enclosed poles. In the former case, also a two
% element vector can be given. The first value is then added to the
% semi-major axis and the second to the semi-minor axis. The convergence is
% optimal if the distance of the background contour to any pole is maximal.
% Please be aware that you have to control the distance to poles outside it
% manually, i.e., when setting the property 'minD' you must check if the
% distances to poles outside the contour are large enough using the plot
% function of RPExpand.
parameters.minD = 1e14;

% Number of integration points for the large contour (background contour).
parameters.nPointsB = 32;

% You can change this value from its default (Inf) to allow for an adaptive
% refinement of the integration points on the contours until the target
% precision is reached. The convergence of the individual contours is
% estimated by comparing the best numerical solution with fewer
% integration points and the sum of all modal components is compared to
% selected direct solutions on the real line if reference solutions are
% given. These are evaluated when expand is called and the property 
% 'referencePoints' (see below) is not empty or when the method 
% 'getReferenceSolution' is called directly. 
parameters.precision = Inf;

% If the property 'referencePoints' is not empty, the scattering problems
% at the specified real frequencies are solved and compared to the
% expansion. 
parameters.referencePoints = [3.8e15 4e15 4.15e15 4.3e15 4.5e15];

% The expansion points define the positions where the expansion is
% evaluated. As additional expansion points do not require solutions to
% additional scattering problems, the number of expansion points does not
% significantly affect the computational costs. 
parameters.expansionPoints = 3.5e15:1e13:4.5e15;

% This parameter is used to define the background contour. The four points
% defined below are used to define a minimal area enclosing ellipse.
parameters.defineBg = [3.5e15 4.5e15]+[-4e14i;2e14i];

% An instance of the class RieszProjection is constructed.
rp = RieszProjection(sc, parameters);

%% Define the contours

% Get the contours. By default, the background contour is an ellipse. For
% more options: execute help rp.getContours in the command window. 
rp.getContours

% Plot contours, expansion points, and poles. If you omit the argument, a
% list of quantities that can be plotted will be printed to the command
% window and you can select with the mouse. 
rp.plot('ComplexPlane')

%% Expand

% Expand the dipole emission. In order to see a list of quantities that 
% can be expanded, omit the argument. Be aware that if you have set the 
% property includeConjugatedPoles to false only quantities are listed which 
% are linear in E. 
plt = {'legendLocation','northwest'};
rp.computeExpansion('DipoleEmission','plot',plt)

%% View results
% You can view the results using the plot function. For a list of
% quantities that can be plotted, execute rp.plot in the command window. 
rp.plot('NormalizedDecayRate')
rp.plot('Error')

%% Expand a cartesian export of the nearfield
plt = {'add2bg',[1 3]}; % add contributions of 1st and 2nd mode to bg
rp.computeExpansion('CartesianExport',rp.referencePoints,'plot',plt)

% In order to remove all files that have been added to project directory
% you can call sc.clean.

%% Adjust plot
% The quantity CartesianExport comes with a custom plot function that
% allows for additional variables. In the same way you can pass variables
% to the method 'plot', using the corresponding parameter of the method
% 'computeExpansion', it is possible to pass arguments to the custom plot
% function, called by 'RieszProjection/plot' using the parameter 'custom'.
% A list of possible options is displayed if you type 
% 'help Scattering/private/viewFields' in the Matlab command window. You
% could for instance have a look at another component: 
rp.plot('CartesianExport','custom',{'component','realy'})

