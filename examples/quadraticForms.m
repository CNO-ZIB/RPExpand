% In this example, expansions of quadratic forms and quantities related 
% to the far field are demonstrated. The former is well suited for this 
% method as no cross terms appear. An expansion of the latter leads to
% finite results as the contour integration allows for circumventing the
% exponential divergence of quasinormal modes outside the resonator.

% Clear the workspace
clearvars
clc

%% Set up the Matlab interface to JCMsuite.
% Activate the third party support and start a daemon. 
% Please refer to the README for details.

% Give the path to the installation directory of JCMsuite and add it to
% the Matlab path. The version should be 4.4.0 or higher.
jcm_root = '~/Software/JCMsuite.5.4.3';
addpath(fullfile(jcm_root, 'ThirdPartySupport', 'Matlab'));

% The directories '@RieszProjection', and '@Scattering' must be accessible.
addpath('../src');

% Start a daemon that handles the jobs submitted to jcmwave_solve.
options = struct('Hostname', 'localhost', ...
                'Multiplicity',8, ...
                'NThreads',2 ... 
                );
            
% Shutdown a possibly running daemon and register a new computer resource.
jcmwave_daemon_shutdown;
jcmwave_daemon_add_workstation(options);

% The working directory for RPExpand can be chosen independently from the
% location of the project files for JCMsuite.
wDir = 'DiamondNanodisk';

%% Construct instance of Scattering and RieszProjection
% Define keys that are passed to jcmwave_solve and call the constructor.
% As in the basic example, the light source is a point source.
keys.finiteElementDegree = 4;
keys.position = [0 145e-9 0];
keys.strength = [1 0 0];

% Path to the project file. The line of code that follows this comment
% assumes, that it is stored in the directory 'scattering' located in the
% working directory 'wDir' defined above. 
projectFile =[wDir filesep 'scattering' filesep 'project.jcmpt'];
sc = Scattering(projectFile, keys, wDir);
rbDir = [wDir filesep 'resultbags_qF'];
if ~exist(rbDir,'dir'), mkdir(rbDir); end
sc.resultbagDir = rbDir; 
rad_keys.gridPointsTheta = 0:5:90;
rad_keys.gridPointsPhi = 0:5:180;
sc.addQuantity('DipoleEmission')
sc.addQuantity('ElectromagneticFieldEnergyFlux')
sc.addQuantity('RadiationPattern2D','RadiationPattern')
sc.addQuantity('RadiationPattern3D','RadiationPattern',rad_keys);
sc.addQuantity('FarFieldIntegral',struct('NA',[0.8 1]),'hidden',true);

%% Construct an instance of RieszProjection
% You can pass a struct with options to the constructor, which will be
% assigned to the corresponding properties. Alternatively, you can set them
% up after construction. For a complete list of public properties, type
% help RieszProjection in the command window.

% In this example, we will expand different quantities including the
% electric field energy and the radiation to the far field. Both are based 
% on functionals quadratic in the electric field. In order to avoid complex 
% conjugation, which does not meet the requirement of holomorphicity for 
% the application of the residue theorem, the complex conjugate of a field 
% E is replaced by the field evaluated at the negative frequency, i.e., 
% E^*(w) -> E(-w), which at real frequencies is the same. The physical 
% quantity is real in the time domain. Hence, E(-w) = E^*(w^*) := E^o(w)
% since E^o(w) has the same poles as E(w) but mirrored at the real line we
% have to consider these 'conjugated poles' as well. With contours being
% symmetric with respect to the real line, no extra scattering problems have
% to be solved to obtain E^o(w).
parameters.includeConjugatedPoles = true;

% All other parameters are the same as in the basic example.
parameters.minD = 1e14;
parameters.nPointsB = 32;
parameters.referencePoints = [3.8e15 4e15 4.15e15 4.3e15 4.5e15];
parameters.expansionPoints = 3.5e15:1e13:4.5e15;
parameters.defineBg = [3.5e15 4.5e15]+[-4e14i;2e14i];

% Construct an instance of the class RieszProjection.
rp = RieszProjection(sc, parameters);

%% Define the contours

% Get and plot the contours.
rp.getContours
rp.plot('ComplexPlane')

%% Expand
% If you execute rp.expand in the command window, a list of quantities is
% printed to the command window which can be expanded with the provided
% post processes. In addition to dipole emission and normalized decay rate,
% which have already been available in the basic example, quantities based 
% on quadratic forms are now included. Since the post process for the
% radiation to the far field includes the option to simultaneously evaluate
% a normalization and since the dipole emission is available as well,
% photon collection efficiency (i.e., the power radiated in a predefined
% numerical aperture (NA) divided by the total power radiated to the upper 
% hemisphere) and the dipole power collection efficiency (i.e., the power 
% radiated in the given NA divided by the total emitted power of the
% source) are appended to the list.

% You can pass all the quantities to be expanded at once. Here, dipole
% emission and electromagnetic field energy flux are essentially the same
% quantities. The first of them is based on the expression
% Gamma(w) = -real(E(w,position).'*strength), where a real-valued strength 
% is defined as the source and therefore fixed even though it is related to 
% E via Maxwell's equations. The second is based on a surface integral 
% of the complex pointing vector S (proportional to cross(E,rot(E)) over a 
% closed surface containing the dipole.
qs = {'NormalizedDecayRate' 'ElectromagneticFieldEnergyFlux' 'Radiation'};
rp.computeExpansion(qs)

%% Plot and control the error
% If you plot the error without specifying the index of a contour, the
% error of the physical quantity is evaluated. As the property 
% 'referencePoints' is not empty, the method 'expand' has called
% 'getReferencesolution' at these points. The reference solutions are used
% to estimate the error. The relative error with respect to them is
% evaluated and the maximum is plotted as a function of the number of
% integration points on the background contour. 
rp.plot('Error')

%% Drop the error below a threshold
% If you want the error to be below a certain threshold, you can set the
% property 'precision' to this value. This property can be a scalar or a 
% vector with two elements. In the latter case, the first value is used as 
% the upper bound of the absolute error and the second one as the upper bound 
% of the relative error. In the former case, both bounds have the same value.
% A call of the method 'expand' will now refine the contours until the 
% precision is reached or the number of integration points on any contour
% exceeds 1024. 
rp.computeExpansion('Radiation','precision',1e-3)

%% Adapt plots
args = {'xunit','nm','legendLocation','northeast','xlim',[416 540]};
rp.plot('NormalizedDecayRate',args{:})

% The unit used for the x-axis is saved and used for other quantities as
% well. If you want to add the unit of a specific quantity, this can be
% done with the keyword argument 'yunit'.
rp.plot('Radiation','yunit','W')

% The default interpreter is LaTeX:
rp.plot('ElectromagneticFieldEnergyFlux','yunit','\frac{J}{s}')

%% Plot combined quantities
rp.plot('PhotonCollectionEfficiency')

%% Visualize radiation pattern
% To get an impression in which direction energy is radiated by a point
% source located in a nanoresonator, you can give discretizations of theta
% and phi in spherical coordinates which are used to create a cartesian 
% grid. The expansion of the radiation pattern is mainly used for 
% visualization purposes, therefore, it is expanded at the reference points 
% instead of the expansion points.
rp.computeExpansion('RadiationPattern2D',rp.referencePoints)

% The plot shows the energy density radiated in the radial direction at the 
% reference frequencies. The plots must be read in spherical coordinates. 
% A large radius means that a lot of energy is radiated in the corresponding 
% direction. The modal contributions can have negative radii. This can be
% understood as the suppression of radiation in the opposite direction. 
% The reference solution is displayed with faces instead of edges. 

%% Radiation pattern in 3D

rp.computeExpansion('RadiationPattern3D',rp.referencePoints)