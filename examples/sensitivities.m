% This example is based on the data publication 
% (https://doi.org/10.5281/zenodo.6614950) and the related paper: 
% 'Computation of eigenfrequency sensitivities using Riesz projections for 
% efficient optimization of nanophotonic resonators' Communications Physics
% volume 5, Article number: 202 (2022) 
% (https://doi.org/10.1038/s42005-022-00977-1).

%% Set up the Matlab interface to JCMsuite
% Activate the third party support and start a daemon. Please refer to the
% README for details.
clearvars; clc

% Provide the path to the installation directory of JCMsuite and add it to
% the Matlab path. The version should be 5.2.0 or higher.
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
wDir = 'NanodiskOnSubstrate';

%% Parameters for JCMsuite and RPExpand

c0 = 299792458; % speed of light

% Parameters for a dielectric nanodisk placed on a three-layer substrate
% published by K. Koshelev, S. Kruk, E. Melik-Gaykazyan, J.-H. Choi, 
% A. Bogdanov, H.-G. Park, and Y. Kivshar, Science 367, 288 (2020).
sample = {'radius',465,'height',635,'alpha',90,'dITO',300,'dSpacer',350};
keys = struct(sample{:});

% Working with the interface to JCMsuite 'Scattering' provided with
% RPExpand the existence of derivatives is determined automatically.
keys.derivativeOrder = 1;
keys.position = [4e-7 3.175e-7 0]; % Position for point evaluation
keys.finiteElementDegree = 5;

% The following parameters adjust the spatial discretization.
keys.hGlobal = 400; % global maximum side length
keys.hLayers = 300; % maximum side length in the layered substrate
keys.hResonator = 200; % maximum side length inside the resonator
keys.hEdges = 2.5; % edge refinement

% Parameters for the initialization of the 'RieszProjection' class. Here,
% we want to define a circular contour around a target frequency.
parameters.defineBg = 2e9*pi*c0/1600; % center of the circle
parameters.minD = parameters.defineBg*1e-2; % radius of the circle

parameters.nPointsB = 8; % number of points on the backround contour
parameters.nMax = 4; % number of expected eigenvalues

pFile = [wDir filesep 'scattering' filesep 'project.jcmpt'];

%% Construct instance of RPExpand and plot contour
sc = Scattering(pFile, keys, wDir);
args = {'position',[4e-7 3.175e-7 0],'component',[0 0 1]};
sc.addQuantity('PointEvaluation',args{:});
rp = RieszProjection(sc, parameters);
rp.getContours('circle');
rp.plot('ComplexPlane')

%% Solve eigenvalue problem
% In the given example only a single eigenvalue is contained in the
% contour. However, there is no restriction, and as long as the eigenvalue
% problem can be solved, the corresponding sensitivities will be added to
% the property 'derivatives'.

ev = rp.computePoles('PointEvaluation');
derivatives = struct2array(rp.derivatives.poles)*1e-10;
ndx = [1 2 4 5 3];
tbl = [1:5;real(derivatives(ndx));imag(derivatives(ndx))];
fprintf(1,'\nParameter %15s %20s\n','Real part','Imaginary part')
fprintf(1,'p%d %22.3f %20.3f\n',tbl)

%% Convergence
% In order to estimate the accuracy we additionally perform a convergence
% study. If you want to recompute the results you have to remove the
% corresponding .mat file.

pMax = 5;
data = cell(1,pMax); 
resultFile = [wDir filesep 'sensitivities.mat'];
if exist(resultFile,'file')
    load(resultFile)
    if length(data)>pMax, data = data(1:pMax); end
end
for it = 1:pMax
    if it<=length(data) && ~isempty(data{it}), continue; end
    keys.finiteElementDegree = it+1;
    sc = Scattering(pFile, keys, wDir);
    rp = RieszProjection(sc, parameters);
    rp.getContours('circle');
    ev = rp.computePoles;
    derivatives = struct2array(rp.derivatives.poles)*1e-10;
    data{it} = [real(derivatives(ndx));imag(derivatives(ndx))];
end
save(resultFile,'data')

%% Plot convergence

ref = data{end};
err_rel = cellfun(@(x){abs((x-ref)./ref)},data(1:end-1));
f2 = figure(2); clf(2); f2.Units = 'centimeters';
f2.Position(3:4) = [9 8];
ax = axes(f2); ax.Box = 'on';
ax.YLabel.String = 'Relative Error (Real Part)'; ax.YLim = [1e-8 5];
ax.XLabel.String = 'd'; ax.YScale = 'log'; ax.NextPlot = 'add';
for it = 1:5
    plot(2:pMax,cellfun(@(x)x(1,it),err_rel),'o-');
    t = sprintf('p_%d',it);
    text(ax,pMax*(1-0.05),8*4^(-it),t,'Color',ax.ColorOrder(it,:))
end
f3 = figure(3); clf(3); f3.Units = 'centimeters';
f3.Position(3:4) = [9 8];
ax = axes(f3); ax.Box = 'on';
ax.YLabel.String = 'Relative Error (Imaginary Part)'; ax.YLim = [1e-8 5];
ax.XLabel.String = 'd'; ax.YScale = 'log'; ax.NextPlot = 'add';
for it = 1:5
    plot(2:pMax,cellfun(@(x)x(2,it),err_rel),'o-');
    t = sprintf('p_%d',it);
    text(ax,pMax*(1-0.05),8*4^(-it),t,'Color',ax.ColorOrder(it,:))
end
