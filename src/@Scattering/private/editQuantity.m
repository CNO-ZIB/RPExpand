function quantity = editQuantity(quantity)
%EDITQUANTITY set defaults and compute final keys if needed

keys = quantity.keys; name = quantity.id;
switch name
    case 'FarFieldIntegral'
        keys = radiationKeys(keys);
        quantity.m = length(keys.NA);
    case 'RadiationPattern'
        keys = radiationPatternKeys(keys);
        quantity.plot = @plotRadiationPattern;
        m = length(keys.gridPointsPhi)*length(keys.gridPointsTheta);
        quantity.m = m;
    case 'referenceDipoleEmission'
        [keys,quantity] = referenceDipoleEmissionKeys(keys,quantity); 
    case 'ModeCoupling'
        [keys,quantity] = modeCouplingKeys(keys,quantity);
    case 'PointEvaluation'
        [keys,quantity] = pointEvaluationKeys(keys,quantity);
    case 'DipoleEmission'
        [keys,quantity] = dipoleEmissionKeys(keys,quantity);
    case 'CartesianExport'
        [keys,quantity] = cartesianExport(keys,quantity);
        quantity.plot = @viewFields;
    case 'ElectromagneticFieldAbsorption'
        keys.omega = [];
end
if ~isfield(keys,'fieldIds'), keys.fieldIds = 1; end
quantity.keys = keys;
end

function keys = radiationKeys(keys)
if ~isfield(keys,'radius'), keys.radius = 1e-5; end
if ~isfield(keys,'NA'), keys.NA = [0.8 1]; end
if ~isfield(keys,'step'), keys.step = 20; end
if ~isfield(keys,'direction'), keys.direction = 'up'; end
% get nodes for Gauss quadrature
gknodes = RieszProjection.nodes15;
nodes = [-gknodes(end:-1:1); 0; gknodes];
% subdevide intervalls in phi
subs = linspace(0,360,ceil(360/keys.step)+1);
subs = [subs(1:end-1); subs(2:end)];
midpt = sum(subs)/2; keys.halfhp = diff(subs)/2;
keys.gridPointsPhi = reshape(nodes*keys.halfhp + midpt,1,[]);
% subdevide intervalls in theta
NA = asind(keys.NA); s = 0; step = keys.step;
if strcmp(keys.direction,'down'), s = 180; NA = 180-NA; end
s1 = linspace(s,NA(1),ceil(abs(s-NA(1))/step)+1);
s2 = s1(end); keys.ndx = length(s1); % end of the first NA
if diff(NA)
    s2 = linspace(NA(1),NA(2),ceil(abs(diff(NA))/step)+1);
end
subs = [s1(1:end-1) s2]; subs = [subs(1:end-1); subs(2:end)];
midpt = sum(subs)/2; keys.halfht = abs(diff(subs)/2);
keys.gridPointsTheta = reshape(nodes*keys.halfht + midpt,1,[]);
end

function keys = radiationPatternKeys(keys)
if ~isfield(keys,'radius'), keys.radius = 1e-5; end
if ~isfield(keys,'fieldIds'), keys.fieldIds = 1; end
if ~isfield(keys,'gridPointsTheta')
    keys.gridPointsTheta = 0:90;
end
if ~isfield(keys,'gridPointsPhi')
    keys.gridPointsPhi = [-90 90];
end
end

function [keys,quantity] = modeCouplingKeys(keys,quantity)
fieldname = '';
for field = fieldnames(keys).'
    c = keys.(field{1});
    if ischar(c) && endsWith(c,'jcmpt')
        fieldname = field{1}; break;
    end
end
if isempty(fieldname)
    error(['Path to a project file (*.jmpt) for propagating'...
           ' mode computations expected']);
end
keys = rmfield(keys,fieldname);

% read project file
if ~exist(c,'file'), error('File %s does not exist',c); end
id = fopen(c, 'r');
jcmpt = fread(id,'*char')';
fclose(id);
mc = quantity.jcmpt;
quantity.jcmpt = c;
[~,fname] = fileparts(c);

s = regexp(jcmpt,Scattering.exportPoints(1:29),'once');
if isempty(s)
    % read mode coupling post process
    msg = sprintf(['\n## post processes for mode overlaps' ...
                   ' with scattering solutions ##\n']);
    jcmpt = [jcmpt msg];
    id = fopen(mc, 'r');
    mc = fread(id,'*char')';
    fclose(id);
    id = fopen(c,'w');
    fwrite(id,[jcmpt sprintf(Scattering.exportPoints) mc],'char');
    fclose(id);
end
if isfield(keys,'extract_polarization')
    p = keys.extract_polarization;
else
    p = [0 0 0];
end
keys.px = p(1);
keys.py = p(2);
keys.pz = p(3);
keys.quantity = 'ElectricFieldStrength';
keys.keepDerivatives = 'no';
keys.addSingularFields = 'yes';
keys.format = 'JCM';
keys.lambda0 = [];
keys.fPath = [fname '_results' filesep 'fieldbag.jcm'];
keys.modeFile = keys.fPath;
keys.outFile = [fname '_results' filesep 'coupling.jcm'];
end

function [keys,quantity] = pointEvaluationKeys(keys,quantity)
if quantity.parent.quantities.DipoleEmission.logical(3)
    error('The quantity dipole emission excludes the point evaluation.')
end
if isfield(keys,'component')
    quantity.component = keys.component;
    keys = rmfield(keys,'component');
elseif ~isfield(quantity,'component')
    quantity.component = [1 1 1];
end
if isfield(keys,'position')
    quantity.position = keys.position;
    keys = rmfield(keys,'position');
elseif ~isfield(quantity,'position')
    error('A position is required for the point evaluation.');
end
if ~isfield(keys,'quantity'), keys.quantity = 'ElectricFieldStrength'; end
keys.keepDerivatives = 'no';
keys.format = 'JCM';
keys.addSingularFields = 'no';
quantity.reference = quantity;
end

function [keys,quantity] = dipoleEmissionKeys(keys,quantity)
if quantity.parent.quantities.PointEvaluation.logical(3)
    error('The quantity point evaluation excludes the dipole emission.')
end
if isfield(keys,'position')
    quantity.position = keys.position;
end
if isfield(keys,'strength')
    quantity.strength = keys.strength;
end
keys = struct();
end

function [keys,quantity] = referenceDipoleEmissionKeys(keys,quantity)
if isfield(keys,'points')
    keys = struct('points',keys.points);
else
    keys = struct('omega',[]);
    quantity.bulkEmission = zeros(0,2);
end
end

function [keys,quantity] = cartesianExport(keys,quantity)
quantity.component = 'realx'; quantity.m = 1;
nm1 = {'keepDerivatives' 'quantity' 'format'};
nm2 = {'NGridPointsX' 'NGridPointsY' 'NGridPointsZ'};
nm3 = {'GridPointsX' 'GridPointsY' 'GridPointsZ'};
vs = {'no' 'ElectricFieldStrength' 'JCM' 100 100 0};
quantity.size = [100 100 1]; fun = @(x,y){sprintf('%s = %d',x,y)};
vs(end-2:end) = cellfun(fun,[nm2(1:2) nm3(end)],vs(end-2:end));
default = cell2struct(vs,[nm1 'X' 'Y' 'Z'],2);
fnames = fieldnames(keys); 
for it = 1:length(fnames)
    switch fnames{it}
        case 'component'
            quantity.component = lower(keys.component);
        case nm1
            default.(fnames{it}) = default.(fnames{it}); 
        case nm2
            n = keys.(fnames{it}); nm = fnames{it}(end);
            default.(nm) = sprintf('%s = %d',fnames{it},n);
            quantity.size(int16(nm)-87) = n;
        case nm3
            n = keys.(fnames{it}); nm = fnames{it}(end);
            default.(nm) = sprintf('%s = %s',fnames{it},mat2str(n));
            quantity.size(int16(nm)-87) = length(n);
        otherwise
            fprintf('Field %s ignored for cartesian export.',fnames{it})
    end
end
quantity.m = prod(quantity.size)*3;
keys = default;

% adapt for versions newer than 6.0.7
if polyval(quantity.parent.version>[5;0;6],2)>4 % 6.0.7 or newer
    keys.localFieldOnly = sprintf('LocalResponseField = %s','yes');
else
    keys.localFieldOnly = sprintf('AddSingularFields = %s','no');
end
end