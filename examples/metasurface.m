% This file is related to the publication 
% "Poles and zeros of electromagnetic quantities in photonic systems" 
% referred to as the paper in later comments.
%#ok<*CTPCT>

%% set up workspace
clearvars % clear workspace 

% provide the path to your intallation of JCMsuite, we used version 5.4.3
jcm_root = '~/Software/JCMsuite.5.4.3';
addpath(fullfile(jcm_root, 'ThirdPartySupport', 'Matlab'));
addpath('../src'); 

% The parameter 'Multiplicity' gives the number of scattering simulations to
% be started in parallel and 'NThreads' defines the number of threads that
% can be used by individual simulations. A single scattering simulation
% requires about 5GB of RAM for the provided accuracy settings.
options = struct('Hostname', 'localhost', ...
                'Multiplicity',8, ...
                'NThreads',2 ... 
                );   
            
% Shut down a possibly running daemon and register a new computer resource.
jcmwave_daemon_shutdown;
jcmwave_daemon_add_workstation(options);

% the project file for scattering simulations
projectFile = ['MetaSurface' filesep 'scattering' filesep 'project.jcmpt'];
% the directory fieldbags are written to
wDir = [pwd filesep 'MetaSurface'];
if ~exist(wDir,'dir'), mkdir(wDir); end
% the directory containing the resultbags
rbDir = [wDir filesep 'resultbags'];
if ~exist(rbDir,'dir'), mkdir(rbDir); end

%% set parameters 
% parameters for JCMsuite project files
keys.finiteElementDegree = 3;
keys.n_glass = 1.5;
keys.A_vec = [1 0 0];
keys.UpperRadius = 250;
keys.Height = 600;

% create interface to JCMsuite
sc = Scattering(projectFile, keys, wDir);
sc.resultbagDir = rbDir; % set directory for resultbags
ft_keys = struct('normalDirection','Z');
ft_settings = {'component',[1 0 0],'diffractionOrder',0};
sc.addQuantity('FourierTransform',ft_keys,ft_settings{:});

% the power flux of the source for normalization
p = 0.5*keys.n_glass*sqrt(RieszProjection.eps0/RieszProjection.mu0);

keys_r = struct('normalDirection','Z'); % keys for reflection
keys_t = struct('normalDirection','-Z'); % keys for transmission

% Reflection and transmission, respectively, must be evaluated from the
% Fourier coefficients in a holomorphic way. The interface to JCMsuite
% comes with the static method 'rt' which does that. In general, if an
% existing post process is used to evaluate a custom quantity, a function
% must be defined, that takes the quantity (a struct with certain fields,
% e.g., 'resultbag') together with tags referring to results in the
% resultbag and returns the quantity of interest. 
r_args = {keys_r,'evaluate',@Scattering.rt,'quadratic',true,'nm',p};
sc.addQuantity('Reflection','FourierTransform',r_args{:});
% The transmission is based on the reflection, i.e., everything is
% identical, except for the keys. 
sc.addQuantity('Transmission','Reflection',keys_t);

% set parameters for RPExpand and create an instance of RieszProjection
parameters.minD = 2.5e14; % radius of the contour
parameters.defineBg = 1.275e15; % center of the contour
parameters.nPointsB = 64; % number of integration points along the contour
parameters.includeConjugatedPoles = true; % include quadratic quantities

 % The property 'nMax' gives the maximum number of poles that can be found
 % inside the contour. It should rather be chosen too large than too small,
 % e.g., if there are 4 poles nMax=8 would lead to good results.
parameters.nMax = 8;

% call the constructor, get contours and plot them
rp = RieszProjection(sc, parameters);
rp.getContours();
y = [-3e14 3.3e14]; x = [9.6e14 15.9e14];
rp.plot('ComplexPlane','contourNumbers',false,'ylim',y,'xlim',x);

%% compute poles and zeros and their derivatives
wn = rp.computePoles('FourierTransform');
zn = rp.computeZeros('FourierTransform');
dw = rp.derivatives.poles;
dz = rp.derivatives.zeros.FourierTransform;

%% update plot
% change legend and colors and add the positions of the zeros (Fig. 2b of
% the paper)
f1 = rp.plot('ComplexPlane');
f2 = figure(2); clf(2)
copyobj(f1.Children,f2);
ax2 = f2.Children(2); ax2.NextPlot = 'add';
ax2.Children(1).Color = 'k';
ax2.Children(1).LineWidth = 2;
ax2.Children(2).MarkerSize = 12;
ax2.Children(2).Color = 'k';
ax2.Children(3).Marker = 'x';
ax2.Children(3).MarkerSize = 6;
ax2.Children(3).LineWidth = 1;
plot(ax2,zn,'bo')
ax2.Children(1).LineWidth = 1;
f2.Children(1).String = {'$\tilde{\omega}_k$'...
    '$\omega_m^{\mathrm{pole}}$' '$\omega_m^{\mathrm{zero}}$'};

%% refine contour for error estimate
rp.nPointsB = 2*rp.nPointsB;
wr = rp.computePoles('FourierTransform');
zr = rp.computeZeros('FourierTransform');
dwr = rp.derivatives.poles;
dzr = rp.derivatives.zeros.FourierTransform;


%% table 1 of the paper
% real and imaginary parts of the pole-zero pair of interest and their
% derivatives with respect to the cone height and the upper radius.
n = 3;
p = [real(wn(n)) imag(wn(n))]; p_r = [real(wr(n)) imag(wr(n))];
err_p = abs(p - p_r)./abs(p_r);
dwH = [real(dw.Height(n)) imag(dw.Height(n))];
dwH_r = [real(dwr.Height(n)) imag(dwr.Height(n))];
err_dpH = abs(dwH - dwH_r)./abs(dwH_r);
dwR = [real(dw.UpperRadius(n)) imag(dw.UpperRadius(n))];
dwR_r = [real(dwr.UpperRadius(n)) imag(dwr.UpperRadius(n))];
err_dpR = abs(dwR - dwR_r)./abs(dwR_r);
z = [real(zn(n)) imag(zn(n))]; z_ref = [real(zr(n)) imag(zr(n))];
err_z = abs(z - z_ref)./abs(z_ref);
dzH = [real(dz.Height(n)) imag(dz.Height(n))];
dzH_r = [real(dzr.Height(n)) imag(dzr.Height(n))];
err_dzH = abs(dzH - dzH_r)./abs(dzH_r);
dzR = [real(dz.UpperRadius(n)) imag(dz.UpperRadius(n))];
dzR_r = [real(dzr.UpperRadius(n)) imag(dzr.UpperRadius(n))];
err_dzR = abs(dzR - dzR_r)./abs(dzR_r);

spec = '%11s%13.4e%12.3e%12.0e%12.0e\n'; spec_ = '%11s%13s%12s%12s%12s\n';
fprintf('\n');
fprintf(spec_,'','Re(u)','Imag(u)','err(Re(u))','err(Im(u))')
fprintf(spec,'pole',[p err_p],'d/dp1',[dwR err_dpR],'d/dp2',[dwH err_dpH]) 
fprintf(spec,'zero',[z err_z],'d/dp1',[dzR err_dzR],'d/dp2',[dzH err_dzH]) 

%% expand the x component of the Fourier coefficient
xx = {11e14 14e14}; yy = {-15e13 5e13}; % define a grid
[x,y] = meshgrid(linspace(xx{:},750),linspace(yy{:},500));
xy = x+1i*y; % the expansion points
% expansion based on the residues
z_m = rp.computeExpansion('FourierTransform',xy(:)); 
% here z(:,:,end) is the background contribution and z{:,:,1:end-1} the
% contributions of the enclosed poles
z = sum(z_m,3);
% z = z(:,:,3); would be the contributions of the third mode only
z = reshape(z(:,1),size(x));

%% plot the phase as a function of complex frequency
% the plot produced with this code corresponds to Fig. 3a in the paper,
% with z = z(:,1,3) only the phase of the third pole can be plotted, which
% would correspond to Fig. 3b
f3 = figure(3); delete(f3.Children); ax3 = axes(f3); ax3.NextPlot = 'add';
args = {'Colormap',hsv,'DisplayRange',[-pi,pi],...
        'Parent',ax3,'XData',x(1,1:end),'YData',y(1:end,1)};
imshow(angle(z),args{:}); ax3.YDir = 'normal';
contour(ax3,x,y,angle(z),20,'Color','k')
plot(ax3,zr,'wo','MarkerSize',6,'LineWidth',1)
plot(ax3,wr,'wx','MarkerSize',6,'LineWidth',1)
plot([1.1 1.4]*1e15,[0 0],'w','LineWidth',1)
ax3.XLim = [min(x(1,1:end)) max(x(1,1:end))];
ax3.YLim = [min(y(1:end,1)) max(y(1:end,1))];
ax3.DataAspectRatio = [1 1 1];
ax3.Visible = 'on'; ax3.TickLabelInterpreter = 'latex';
ax3.FontSize = 14;
cb = colorbar(ax3);
cb.TickLabelInterpreter = 'latex';
cb.Ticks = [-pi -pi/2 0 pi/2 pi];
cb.Label.Interpreter = 'latex'; cb.Label.String = 'Arg($q$)';
tLabels = {'$-\pi$' '$-\frac{\pi}{2}$' '$0$' '$\frac{\pi}{2}$' '$\pi$'};
cb.TickLabels = tLabels; cb.FontSize = 14;
ylabel(ax3,'Im($\omega$) [Hz]','Interpreter','latex')
xlabel(ax3,'Re($\omega$) [Hz]','Interpreter','latex')
nm = {'$\omega_m^{\mathrm{pole}}$' '$\omega_m^{\mathrm{zero}}$'};
leg = legend(ax3.Children([2 3]),nm,'Interpreter','latex');
leg.Location = 'northwest';
leg.Color = [0.7 0.7 0.7]; leg.ItemTokenSize = [15;9];

%% expand physical observalbes
% Since reference points are defined, the reference solutions is evaluated
% automatically
w0 = linspace(1.1e15,1.45e15,1001).'; % define frequency range or interest
rp.expansionPoints = w0;
rp.referencePoints = w0(1:100:end);
reflection = rp.computeExpansion('Reflection','plot',true);
transmission = rp.computeExpansion('Transmission','plot',true);

%% check energy conservation
f_e = figure(123); clf(123); ax_e = axes(f_e); ax_e.NextPlot = 'add';
r = sum(real(reflection),3);
t = sum(real(transmission),3);
s = r+t;
plot(ax_e,rp.expansionPoints,r,'DisplayName','reflection')
plot(ax_e,rp.expansionPoints,t,'DisplayName','transmission')
plot(ax_e,rp.expansionPoints,s,'DisplayName','sum')
ax_e.YLim = [-0.1 1.3];
legend(ax_e)
