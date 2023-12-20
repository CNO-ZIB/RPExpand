classdef Scattering < handle
    %SCATTERING interface to FEM solver JCMsuite
    %   This class parallelizes the evaluation of scattering problems with
    %   all parameters fixed except the frequency. An instance sc of this
    %   class can be used as a callable which takes a cell array of 
    %   doubles and returns tags representing results of scattering 
    %   problems. Based on these results the physical quantities to be
    %   expanded are evaluated in corresponding postprocesses.
    %   If the property 'keys' has a field 'position', it is assumed to
    %   be the position of a dipole located within the computational
    %   domain, and the field value of the scattered field at its position
    %   is exported. Please refer to the <a href="matlab:web( ...
    %   ['https://docs.jcmwave.com/JCMsuite/html/ParameterReference' ...
    %   '/d77be858da38e480d80e9f59c1a3f25f.html'],'-browser')"
    %   >parameter reference</a> 
    %   of JCMsuite if you want to learn more about the definition of  
    %   scattering problems. The project file, e.g., project.jcmpt, must
    %   contain the section PML { %(pml)s } and the primitive 
    %   FiniteElementDegree = %(finiteElementDegree)e for parameter
    %   substitution. Please refer to the documentation of the
    %   <a href="matlab:web(['https://docs.jcmwave.com/JCMsuite/html/' ...
    %   'MatlabInterface/index.html'],'-browser')" >matlab interface</a>.
    %   Your custom parameters for keyword substitution can be passed via
    %   the property 'keys' of this class which must at least have the
    %   field 'finiteElementDegree'. Apart from the project file, you need 
    %   the files <a href="matlab:web(['https://docs.jcmwave.com/' ...
    %   'JCMsuite/html/ParameterReference/b61236968b3822be5ffbfee' ...
    %   '6564f23da.html'],'-browser')"
    %   >layout.jcm</a> or layout.jcmt, <a href="matlab:web(['https:' ...
    %   '//docs.jcmwave.com/JCMsuite/html/ParameterReference/' ...
    %   '3df274a2924c89630ff2393cc22b686e.html'],'-browser')"
    %   >materials.jcm</a> or materials.jcmt
    %   and <a href="matlab:web(...
    %   ['https://docs.jcmwave.com/JCMsuite/html/ParameterReference/' ...
    %   'f3e666a5067147d3cd45b67773bb77ae.html'], '-browser')"
    %   >sources.jcmt</a>. The latter has to contain the 
    %   primitive Omega = %(omega)e and, in case of a <a href="matlab:
    %   web(['https://docs.jcmwave.com/JCMsuite/html/Parameter' ...
    %   'Reference/67f12453f439cffbc214bb0a2362a853.html'],'-browser')"
    %   >point source</a>,
    %   additionally Position = %(position)e and Strength = %(strength)e
    %   are required which then must also be fields of the property 'keys'.
    %
    %SCATTERING Properties:
    %   keys - Struct passed to JCMsolve for parameter substitution
    %   These keys are used for solving scattering problems. Keys for post
    %   processes must be provided with the quantities.
    %   projectFile - Path to a JCMsuite template file
    %   The file must end with .jcmpt and its location must be a project
    %   directory containing all files needed for the scattering problem.
    %   derivativeNames - Names of derivative parameters 
    %   If the derivative order in the file project.jcmp is set to one,
    %   as soon as results with derivative information are processed, the
    %   names of the derivative parameters are stored.
    %   workingDir - The directory in which results will be saved
    %   eigenvalues - Eigenvalues are loaded from the 'evsDir'
    %   pml - File containing the pml settings
    %   resultbagDir - Directory for resultbags
    %   evsDir - Directory containing results of a resonance problem
    %
    %JCMSCATTERING Methods: 
    %   Scattering - Constructor
    %   clean - Delete Files that have been created using the class
    %   'RieszProjection'
    %   dipoleBulkEmission - Get the dipole emission in bulk material
    %   addPostProcess - Define quantity based on post process
    %   customizeQuantity - Redefine existing quantity
    %
    %DEPENDENCIES:
    %   The JCMwave third party support has to be set up for this class.
    %   And you have to run a daemon in order to apply post processes.
    %   Please refer to the documentation of the <a href="matlab:
    %   web(['https://docs.jcmwave.'com/JCMsuite/html/MatlabInterface' ...
    %   '/index.html'],'-browser')">matlab interface</a> to set up
    %   the third party support and to the homepage of <a href="matlab:
    %   web('https://jcmwave.com/jcmsuite/jcmsuite-documentation', ...
    %   '-browser')">JCMwave</a> for
    %   installation instructions.
    %
    %   See also RIESZPROJECTION

% This file is part of the software RPExpand
% Copyright: 2021 Zuse Institute Berlin
% Authors: Fridtjof Betz, Felix Binkowski
% Updated: April-2023
    
    
    properties (SetAccess=private)
        keys (1,1) struct {mustHaveField} = ...
            struct('finiteElementDegree',0);
        
        projectFile (1,:) char {mustBeValidFile};
        
        derivativeNames (1,:) cell;
        
        workingDir (1,:) char {mustBeValidDirectory};
        
        eigenvalues (:,1) double;
    end
    
    properties
        pml (1,:) char;
        
        resultbagDir (1,:) char {mustBeValidDirectory};
    end
    
    properties (Dependent)
        evsDir (1,:) char {mustBeValidDirectory};
    end
    
    properties (Hidden)
        % projectDir - The directory where the project files are located
        projectDir (1,:) char {mustBeValidDirectory};
        
        % postproDir - The directory where post processes are located
        postproDir (1,:) char {mustBeValidDirectory};
        
        % scatteringDir - Working directory for the scattering problems
        scatteringDir (1,:) char {mustBeValidDirectory}; 
        
        % referenceDir - Working directory for the reference solutions
        referenceDir (1,:) char {mustBeValidDirectory};
        
        % resultbags - Resultbags of the scattering problems
        % resultbags{1}: contour, resultbags{2}: reference points
        resultbags (1,2) cell;
        
        % bulkEmission - Dipole bulk emission of the given setting
        bulkEmission (:,2) double; 

        % diffractionOrder - diffraction order (Fourier transform)
        diffractionOrder = 0;
        
        % pprbs - resultbags for some postprocesses
        % used by the method RieszProjection/viewFields
        pprbs = struct;
    end
    
    properties (Access=private)
        % pointSource - If true, the source file contains a single dipole
        pointSource (1,1) logical = false;
        
        % see dependent properties
        evsDir_
        
        % quantities - based on postprocesses
        quantities (1,1) struct;
        
        % resonance - contains data required for resonance expansions
        % resonance(1:3): position of the dipole
        % resonance(4:6): strength of the dipole
        % resonance(7:2:end-1): field exports
        % resonance(8:2:end): eigenfrequencies
        resonance (:,1) double;
        tags (:,1) cell;
        % derivatives - true if derivatives are available
        derivatives = false;
        
        % version of JCMsuite
        version double;
        
        % the instance of the class 'RieszProjection'
        parent RieszProjection = []; 
    end
    
    methods
        function sc = Scattering(projectFile, keys, wDir)
            %SCATTERING Construct an instance of this class
            %   SCATTERING(projectFile, keys) constructs an object with
            %   working directories located in the parent directory of the
            %   .jcmpt file 'projectFile', which must define a scattering
            %   problem. The struct keys contain fields that are
            %   used for the parameter substitution of JCMsuite. It must a
            %   least contain the field 'finiteElementDegree'.
            %   SCATTERING(projectFile, keys, wDir) set the working
            %   directory to wDir. The working directory contains all files
            %   and directories created and written to by this object.
            
            if nargin<3, wDir = ''; end
            % make sure all keys exist, that could be relevant
            fields = Scattering.fields;
            for it = 1:length(fields)
                if isfield(keys,fields{it})
                    error(['The fieldname %s is reserved for ', ...
                        'internal use.'], fields{it})
                else
                    keys.(fields{it}) = [];
                end
            end
            sc.keys = keys;
            sc.setPaths(projectFile,wDir);
            pD = dir(sc.projectDir);
            % check if a single point source is contained
            sF = pD(startsWith({pD.name},'sources.jcm'));
            s = fileread([sF(1).folder filesep sF(1).name]);
            c = count(s,'Source {')==2 && count(s,'PointSource {')==1;
            sc.pointSource = c;
            addDefaultPostProcesses(sc);
            
            global JCMwaveGlobal;
            if isfield(JCMwaveGlobal,'version')
                sc.version = str2double(split(JCMwaveGlobal.version,'.'));
            end
        end
        
        function clean(sc,deleteResultbags)
            %CLEAN Delete files
            %   CLEAN(sc) Removes all scattering problems in the working
            %   directories integration_points and reference_points as well
            %   as the directory 'superpositions' and 'bulk_emission' if
            %   they exist. The file 'pml.log' is deleted too. If you want
            %   to use a custom file that you do not want to be deleted,
            %   you should place it outside the project directory and set
            %   the property 'pml' to the file path.
            %
            %   CLEAN(sc, deleteResultbags) If deleteResultbags is false
            %   the directory containing the resultbags is kept. The
            %   default is true.
            
            if nargin<2, deleteResultbags = true; end
            if exist(sc.scatteringDir,'dir')
                rmdir(sc.scatteringDir,'s');
                mkdir(sc.scatteringDir);
            end
            if exist(sc.referenceDir,'dir')
                rmdir(sc.referenceDir,'s');
                mkdir(sc.referenceDir);
            end
            if deleteResultbags && exist(sc.resultbagDir,'dir')
                rmdir(sc.resultbagDir,'s');
                mkdir(sc.resultbagDir);
                sc.resultbags = cell(1,2);
                sc.pprbs = struct;
            end
            sp = [sc.workingDir filesep 'superpositions'];
            if exist(sp,'dir'), rmdir(sp,'s'); end
            bk = [sc.workingDir filesep 'bulk_emission'];
            if exist(bk,'dir'), rmdir(bk,'s'); end
            if exist(sc.pml,'file') && contains(sc.pml,sc.workingDir)
                delete(sc.pml); 
            end
        end
        
        emission = dipoleBulkEmission(sc,omega,bE);
        
        function set.bulkEmission(sc,bE)
            if size(sc.bulkEmission,1)>size(bE,1)
                idx = ismembertol(bE(:,1),sc.bulkEmission(:,1));
                bE = sortrows([bE(~idx,:); sc.bulkEmission]);
            end
            sc.bulkEmission = bE;
        end
        
        function evsDir = get.evsDir(sc)
            evsDir = sc.evsDir_;
        end
            
        % add quantity based on a given post process
        addPostProcess(sc,name,varargin)
        
        % adapt existing quantities 
        customizeQuantity(sc,name,varargin)
        
        viewEigenmodes(sc,varargin)
        
        function set.evsDir(sc,evsDir)
            pF = java.io.File(evsDir);
            if ~pF.isAbsolute
                pF = java.io.File([pwd filesep evsDir]);
            end
            pF = pF.getCanonicalFile;
            sc.evsDir_ = char(pF.getPath);
            evs = jcmwave_load([evsDir filesep 'eigenvalues.jcm']);
            sc.eigenvalues = evs.eigenmode;
        end
    end
    
    methods (Hidden)
        function varargout = subsref(sc, S)
            %SUBSREF Subscripted reference
            %   To allow for calls sc(contours) of an instance sc of this 
            %   class the subsref method has to be redefined. sc(contours)
            %   calls the method f which can be redefined and is supposed
            %   to solve scattering problems on given contours in the
            %   frequency plane.
            %
            %   See also Scattering.F, SUBSREF
            if strcmp(S(1).type, '()') && iscell(S(1).subs{1})
                varargout{1} = sc.f(S(1).subs{:});
            else
                if nargout>0
                    varargout{1:nargout} = builtin('subsref', sc, S);
                else
                    builtin('subsref', sc, S);
                end
            end
        end
        
        % This function is called by subsref and can be redefined for 
        % different scattering problems. The input parameter contours are 
        % expected to be a cell of size=(1,:) containing doubles of 
        % shape=(1,:).
        out = f(sc, varargin)
        
        out = sum(sc, tags, weights, omega0, outName)
        
        [ids,tags] = submit(sc, tags, varargin)
        
        v = collect(sc, ids, tags, q, last)
        
        varargout = export(sc, tags, pointListOrKeys, quantity)
        
        out = contributions(sc, q, w, w_n, v_n, keys)
    end
    
    methods (Access = ?RieszProjection)
        quantitiesForExpansion(sc,rp)
    end
    
    methods (Static)
        varargout = load(filename)
        out = rt(quantity,tags)
    end
    
    methods (Access=private)
        function sc = setPaths(sc, projectFile, wDir)
            %SETPATHS Create directories to work in
            %   SETPATHS(sc, projectFile) Creates directories located in
            %   the parent directory of projectFile.
            
            a = filesep;
            pF = java.io.File(projectFile);
            if ~pF.isAbsolute
                pF = java.io.File([pwd a projectFile]);
            end
            pF = pF.getCanonicalFile;
            sc.projectFile = char(pF.getPath);
            pD = char(pF.getParent);
            sc.projectDir = pD;
            if isempty(wDir), wDir = pD; end
            wF = java.io.File(wDir);
            if ~wF.isAbsolute, wF = java.io.File([pwd a wDir]); end
            wDir = char(wF.getCanonicalPath);
            sc.workingDir = wDir; 
            scDir = [wDir a 'integration_points'];
            if ~exist(scDir,'dir'), mkdir(scDir); end
            refDir = [wDir a 'reference_points'];
            if ~exist(refDir,'dir'), mkdir(refDir); end
            resDir = [wDir a 'resultbags'];
            if ~exist(resDir,'dir'), mkdir(resDir); end
            sc.scatteringDir = scDir;
            sc.referenceDir = refDir;
            sc.resultbagDir = resDir;
            pD = mfilename('fullpath'); n = strfind(pD, a);
            pD = pD(1:n(end-2)-1);
            sc.postproDir = [pD a 'postprocesses'];
            if ~exist(sc.postproDir,'dir')
                warning(['The directory containing the post processes '...
                    'is expected to be found in the directory '...
                    'containing the folder with the definition of this '...
                    'class. You must either copy it there or set the '...
                    'hidden property postproDir of this class.'])
            end
        end
        
        function varargout = tag2path(sc, tags)
            %TAG2PATHS convert tags into file paths
            
            files = {}; res = [];
            if exist(tags{1},'file')
                files = tags; 
            else
                isref = ~isempty(sc.resultbags{2}) && ...
                    isfield(sc.resultbags{2}.results_,tags{1});
                res = sc.resultbags{1+isref}.results_;
                if isfield(res,tags{1})
                    files = cellfun(@(x){res.(x).result{1}.file},tags);
                end
            end
            varargout{1} = files;
            if nargout>1
                w = [];
                if ~isempty(res)
                    w = cellfun(@(x)res.(x).keys.omega,tags);
                end
                varargout{2} = w;
            end 
        end
    end
    
    properties (Constant)
        eps0 = 8.85418781762039e-12; % permittivity in vacuum
        mu0 = 4e-7*pi; % permeability in vacuum
        c0 = 299792458; % speed of light
    end
    
    properties (Constant, Access=private)
        bags = {'integration_points.mat' 'reference_points.mat'};
        fields = {'addSingularFields' 'px' 'py' 'pz' 'quantity' 'omega'};
        
        pP = [ ... dipole bulk emission post process
            'PostProcess {\n' ...
            '  BulkEmission {\n' ...
            '    Omegas=%%(omega)e\n' ...
            '    ProjectPath="%s"\n' ...
            '    OutputFileName="bulk_emissions.jcm"\n' ...
            '  }\n' ...
            '}\n'...
            ];
        
        superposition = [ ... superposition post rprocess
            'PostProcess { \n' ...
            '  Superposition { \n' ...
            '    OutputFileName = "%%(oPath)s"\n' ...
            '%%(fieldbags)s\n' ...
            '    Omega = %%(omega)e\n' ...
            '  } \n' ...
            '}\n'];

        fieldBag = [ ... field bags for superposition post process
            '    FieldBag { \n' ...
            '      FieldBagPath = "%s"\n' ...
            '      Weights = %s \n' ...
            '    }\n'];
        
        exportPoints = [ ... point list
            'PostProcess {\n' ...
            '  ExportFields {\n' ...
            '    FieldBagPath = "%%(fPath)s"\n' ...
            '    OutputFileName = "%%(oPath)s"\n' ...
            '    OutputQuantity = "%%(quantity)s"\n' ...
            '    %%(localFieldOnly)s\n' ...
            '    KeepDerivatives = %%(keepDerivatives)s\n' ...
            '    Format = %%(format)s\n' ...
            '    Cartesian {\n' ...
            '      GridPointsX = %%(px)e\n' ...
            '      GridPointsY = %%(py)e\n' ...
            '      GridPointsZ = %%(pz)e\n' ...
            '    }\n' ...
            '  }\n' ...
            '}\n'];
                    
        exportCartesian = [ ... cartesian
            'PostProcess {\n' ...
            '  ExportFields {\n' ...
            '    FieldBagPath = "%%(fPath)s"\n' ...
            '    OutputFileName = "%%(oPath)s"\n' ...
            '    OutputQuantity = "%%(quantity)s"\n' ...
            '    Cartesian {\n' ...
            '      %s\n' ...
            '      %s\n' ...
            '      %s\n' ...
            '    }\n' ...
            '  }\n' ...
            '}\n'];
        
        normalizationConstant = [ ... QNMs normalization
            'PostProcess {\n' ...
            '  DensityIntegration {\n' ...
            '    OutputFileName = "out1.jcm"\n' ...
            '    OutputQuantity = HolomorphicElectricFieldEnergy\n' ...
            '    FieldBagPath = "%%(fPath)s"\n' ...
            '  }\n' ...
            '}\n' ...
            'PostProcess {\n' ...
            '  DensityIntegration {\n' ...
            '    OutputFileName = "out2.jcm"\n' ...
            '    OutputQuantity = HolomorphicMagneticFieldEnergy\n' ...
            '    FieldBagPath = "%%(fPath)s"\n' ...
            '  }\n' ...
            '}\n'];
        
        source = [ ... dipole bulk emission post process
            'SourceBag {\n' ...
            '  Source {\n' ...
            '    ElectricCurrentDensity {\n' ...
            '      PointSource {\n' ...
            '        Omegas=%%(%s)e\n' ...
            '        Position=%%(position)e\n' ...
            '        Strength=%%(strength)e\n' ...
            '      }\n' ...
            '    }\n' ...
            '  }\n' ...
            '}\n'...
            ];   
    end
end

function mustBeValidFile(filepath)
if ~exist(filepath, 'file') && ~isempty(filepath)
    error('The file \"%s\" does not exist', filepath);
end
end

function mustBeValidDirectory(filepath)
if ~exist(filepath, 'dir') && ~isempty(filepath)
    error('The directory \"%s\" does not exist', filepath);
end
end


function mustHaveField(keys)
if ~isfield(keys,'finiteElementDegree')
    error('The struct keys must contain the field finiteElementDegree');
end
end
