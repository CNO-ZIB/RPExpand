% This example is based on the data publication 
% (https://doi.org/10.5281/zenodo.7376555) and the related paper: 
% 'Resonance expansion of quadratic quantities with regularized quasinormal
% modes' Physica Status Solidi A, volume 220, issue 7 (2023)
% (https://doi.org/10.1002/pssa.202370013).

%% Set up the Matlab interface to JCMsuite
% Activate the third party support and start a daemon. Please refer to the
% README for details.
clearvars; clc

% Provide the path to the installation directory of JCMsuite and add it to
% the Matlab path. The version should be 5.2.1 or higher. Here, for
% versions newer than 6.0.0 the mode coupling post process does not work as
% expected.
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
            
% Shutdown a possibly running daemon and register a new computer resource.
jcmwave_daemon_shutdown;
jcmwave_daemon_add_workstation(options);

% The working directory for RPExpand can be chosen independently from the
% location of the project files for JCMsuite.
wDir = 'CircularBraggGrating';

%% Prepare the parameters

skeys.finiteElementDegree = 4; % finit element degree

% Material constants
skeys.n_CBG = 3.39;
skeys.n_HSQ = 1.4;
skeys.n_ITO = 1.25;
skeys.k_ITO = 0.07;
skeys.n_SiO2 = 1.45;
skeys.n_Au = 0.4;
skeys.k_Au = 8.95;
skeys.n_Core = 1.4885;
skeys.n_Cladding = 1.4468;

% Geometry parameters
skeys.nRings = 10; % number of rings
skeys.t_CBG = 255; % height of the grating
skeys.t_SiO2 = 400; % heigth of the SiO2 layer
skeys.t_HSQ = 455; % height of the HSQ layer
skeys.t_ITO = 50; % height of ITO layer
skeys.fiberLength = 500; % length of fiber above the ITO layer
skeys.R = 320; % radius of the inner disk
skeys.W = 145; % gap width
skeys.w = 335; % ridge width
skeys.displacement = -(skeys.fiberLength+skeys.t_ITO+skeys.t_HSQ);
skeys.width_total = skeys.R+skeys.nRings*(skeys.w+skeys.W)+1000; 
skeys.r_Core = 900; % UHNA3

% FEM mesh constraints
skeys.h_global = 200;
skeys.h_center = 100;
skeys.h_flux = 25;
skeys.h_corners_center = 5;
skeys.h_rings = 100;
skeys.h_corners = 10;

% Source
skeys.strength = [2e-7 0 0]; 
skeys.position = [0 (skeys.displacement+skeys.t_CBG/2)*1e-9 0];

scattering_jcm =[wDir filesep 'scattering' filesep 'cbg.jcmpt'];

% Resonance problem
rkeys = rmfield(skeys,{'strength','position'});
rkeys.guess = round(2*pi*299792458/1310e-9,4,'significant');
rkeys.nEvs = 30;
rkeys.cylinder = 'yes'; % get linearly polarized eigenmodes

%% Parameters for RPExpand

% We are interested in quantities depending quadratically on the electric
% field strength E.
parameters.includeConjugatedPoles = true;

% Eigenvalues and normalized eigenvectors are used for the expansions.
parameters.useEigenvectors = true;

% Define the frequency range of interest and some reference points.
parameters.expansionPoints = linspace(1.313e15,1.467e15,1541);
parameters.referencePoints = parameters.expansionPoints(1:70:end);

% Instead of using a background contour, for this example, we interpolate
% the background contributions based on a few scattering problems at real
% frequencies.
parameters.nInterpolationNodes = 8;

% The quantities to be expanded:
qs = {'DipoleEmission','AbsorptionITO','CouplingEfficiency'};

%% Solve eigenvalue problem
% We use the Arnoldi method, implemented in JCMsuite, as it provides both
% the normalized eigenvectors and the eigenvalues. In principle, it is
% possible to get both with the Riesz projection method too but the
% normalization is not yet implemented.

evDir = [wDir filesep 'resonance'];
id = jcmwave_solve([evDir filesep 'cbg.jcmp'], rkeys);
jcmwave_daemon_wait(id);
copyfile([wDir filesep 'resonance' filesep 'pml.log'],wDir);

%% 
% Get interface to JCMsuite
sc = Scattering(scattering_jcm, skeys, wDir);
sc.evsDir = [evDir filesep 'cbg_results'];

% Absorption only takes place in the ITO layer, which is associated with the
% DomainId four. Furthermore, we are interested in the coupling to the
% fiber above the ITO layer.
fkeys.n_Core = skeys.n_Core; % refrective index of the core
fkeys.n_Cladding = skeys.n_Cladding; % refrective index of the cladding
fkeys.finiteElementDegree = skeys.finiteElementDegree;
fkeys.r_Core = skeys.r_Core; % radius of the core

% We are interested in how the scattering solution couples to the propagating
% modes of a single mode fiber. Its modes are evaluated based on the
% project files in the directory 'UHNA'.
fkeys.modeFile = [pwd filesep 'UHNA' filesep 'fiber.jcmpt'];
fkeys.position = [0 0 -400e-9]; % position where the overlap is evaluated
akeys = struct('domainIds',4); % the absorption is evaluated in domain 4
sc.addQuantity(qs{2},'ElectromagneticFieldAbsorption',akeys);
sc.addQuantity('FiberCoupling','ModeCoupling',fkeys);
sc.addQuantity('DipoleEmission')

% Here we define the coupling efficiency which is the ratio between the
% energy coupled to the fiber and the total energy emitted by the dipole.
eta.quadratic = true;
eta.parents = {'FiberCoupling' 'DipoleEmission'};
eta.evaluate = @(cp,dpe)cp(:,1,:)./real(sum(dpe,3));

% Construct an instance of the class RieszProjection.
rp = RieszProjection(sc, parameters);
rp.quantities.CouplingEfficiency = eta;

selection = abs(sum(rp.expansionPoints([1 end]))/2-rp.poles)<1e14;
rp.selectedPoles = num2cell(find(selection));

%% Plot contours
args = {[-8e13,8e13],'xlim',[1.25e15 1.65e15]};
f1 = rp.plot('ComplexPlane','ylim',args{:});

%% Expand
% Expand quantities

% For the plot contributions of these modes are added to the background.
add2bg = [1 3 4 5 7 10]; xl = [1.313e15 1.467e15]; % x axis limit
plt = {'legendLocation','northwest','add2bg',add2bg,'xlim',xl};
rp.computeExpansion(qs,'plot',plt)
