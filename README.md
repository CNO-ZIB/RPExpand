# RPExpand

The software RPExpand is designed for efficient evaluations of
electromagnetic quantities based on rational approximants whose summands
can be attributed to eigenmodes of the physical system. Based on contour 
integration methods [1], eigenfrequencies are determined that, together 
with corresponding residues, are used to reconstruct the target quantity
at any frequency inside the contour [2]. Furthermore, its zeros can be 
retrieved and, optionally, derivative information along the contour can be
used to compute derivatives of eigenfrequencies and zeros [3,4]. 

The target quantity does not need to be linear in the electric or magnetic
field but can be based on quadratic forms and derived from the far
field [5]. The results can be visualized with a built-in plot function. For
a comprehensive control of the error, expansions are compared to direct 
solutions at selected frequencies and to results with reduced numbers of 
integration points.

The number of modes present in the spectral region of interest and
coupling to the specified source can be large. Then, it can be preferable to use
eigenvalues and corresponding normalized eigenvectors from external 
solvers to perform expansions based on quasinormal modes (QNMs) and 
polynomial interpolation [6]. The resulting modal contributions are 
equivalent to the residue based approach described above and quantities
quadratic in the electric field can be expanded directly and without cross 
terms using a regularization method [7]. The regularization is based on
scattering simulation at the complex conjugated eigenfrequencies and 
provides accurate results at any distance from the resonator.

One advantage of contour integral methods is parallelizability. This has
been extensively exploited in the design of the interface to the finite 
element method (FEM) solver [JCMsuite](https://jcmwave.com/). The interface 
handles the parallel submission of jobs keeping the interaction with the 
solver as easy as possible while allowing advanced users to access its
full range of functionalities. 

## Installation

Get a copy of the two class folders and the directory containing post
processes of JCMsuite used to derive the target quantities from the 
electric field:

  - @RieszProjection 
  - @Scattering (interface to JCMsuite)
  - postprocesses (required by the interface to JCMsuite)

and add their parent directory to your MATLAB path or copy them to the
directory containing your scripts.

### Requirements

RPExpand has been tested with MATLAB R2018b under Linux but should be
independent of the platform. 

A custom interface must have the form of a callable object (e.g., a
function handle). A rather complex example is the class 'Scattering'
which is the interface to JCMsuite whose version should be 4.4.0 or higher.
A much simpler example is provided in the file 'quantumExample.m'. There, 
the interface is a pure MATLAB function. Input and output of such a custom 
interface are described in more detail in the section 'Notes on a usage 
independent of JCMsuite' at the end of this README.

## Documentation

Comments explaining the usage of classes and functions can be displayed
using the help function. E.g., `help RieszProjection` will print a short
description and a list of properties and methods, which are linked to more
detailed instructions on the corresponding items. 

Detailed examples are provided to make you familiar with the different 
problem types and can serve as starting points for setting up your own 
projects.

Additionally, subsequent sections of this README describe the requirements 
for the interface to the FEM solver JCMsuite and list predefined 
quantities. As the class 'RieszProjection' is not restricted to be used in 
combination with JCMsuite the last section hints on its use independent of 
the interface 'Scattering'. Indeed, the only requirement are functions, 
that can be analytically continued to the complex plane. In the context of 
optics, these functions are optical observables based on solutions of 
Maxwell's equations, solved at complex frequencies. A pure MATLAB example, 
which is not restricted to be used in combination with JCMsuite, is 
presented in the script 'quantumExample.m'.

## Concept

The first step is always the definition of a contour or a selection of poles
to be considered. This is done with the method 'getContours' and different 
properties of the class 'RieszProjection'. The corresponding setup can be
reviewed calling `rp.plot('ComplexPlane')`. Once this figure exists, it will 
be kept up to date if you change properties that effect the displayed quantities.

The function values of interest at the integration points are not evaluated 
unless you call a method starting with 'compute', e.g., 'computeExpansion' or
'computePoles'. An exception is the method 'plot'. If you request the plot of a
quantity that has not yet been expanded, the method 'computeExpansion' will be
called internally. All the visible data, i.e., data that is shown in a figure,
will be kept up to date. 

The expansions are not stored in a property. If you call `rp.computeExpansion(...)`
without a return value, the result will be displayed in a figure. Alternatively
you can allways get the corresponding data as a return value, e.g.,
`value = rp.computeExpansion(...)`. Once the values along the contour have been
computed, subsequent calls of the method 'computeExpansion' for the corresponding
quantity require negligible computational effort. The returned 
expansions have the size (n,d,m) where n is the number of expansion points, 
d a custom data dimension (usually 1) and m the number of contributions. The 
function values along the expansion points can be reconstructed as the sum 
`q(w0) = sum(expansion,3)`.

## The interface to JCMsuite

The current release can be downloaded following this 
[link](https://installation.jcmwave.com/).
Please follow the installation and activation instructions provided there.
A 14-day trial license is available free of charge. It does 
not allow for the parallel computation of different scattering problems, 
however, multithreading is supported. Nevertheless, you need to start a 
daemon.

### Start a daemon

As a prerequisite for the usage of RPExpand the MATLAB interface to
JCMsuite must be set up and a daemon has to be started. A minimal 
example given a personal laptop or desktop machine is:

```matlab
% Set up thirdparty support
jcm_root = '/path/to/local/installation/JCMsuite.x.x.x'
addpath(fullfile(jcm_root,'ThirdPartySupport','Matlab'))

% shutdown a possibly running daemon
jcmwave_daemon_shutdown();

% register a new resource
options = struct(...
                 'Hostname', 'localhost', ...
                 'Multiplicity', 1, ...
                 'NThreads', 2, ...
                 );
jcmwave_daemon_add_workstation(options);
```

If you have access to a remote machine via ssh, you can provide the 
corresponding hostname, user name, etc. For information about additional 
parameters and the use of a queue or cluster, please refer to the
[daemon command reference](https://www.docs.jcmwave.com/JCMsuite/html/MatlabInterface/6b7e17ce176f4ffde673cced4f12eeee.html).

If you have access to a license that allows for the parallel submission of
different jobs, you can increase the multiplicity to a value larger than 
one depending on the available kernels of the computer. The parameter 
'NThreads' refers to the number of threads per job.

### Custom project requirements

This section gives a brief overview of the creation of a directory hosting 
a scattering project of JCMsuite compatible with the present code. It must
contain the following files: 
- [project.jcmpt](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/db9062933554f66e7fb21c46d53fcca2.html)
- [sources.jcmt](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/f3e666a5067147d3cd45b67773bb77ae.html)
- [materials.jcm(t)](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/3df274a2924c89630ff2393cc22b686e.html)
- [layout.jcm(t)](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/b61236968b3822be5ffbfee6564f23da.html)

The trailing 't' at the file extensions marks template files, which are 
used for parameter substitution. Whereas in the last two cases, it is
optional, the first two files must contain predefined parameters. A minimal 
project file is: 

```
Project {
    Electromagnetics {
        TimeHarmonic {
            Scattering {
                PML {
                    %(pml)s
                }
                FieldComponents = Electric
                Accuracy {
                    FiniteElementDegree = %(finiteElementDegree)e
                }
            }
        }
    }
}
```
The parameters 'pml' and 'finiteElementDegree' will be set by the program. 
The former, as all integration points should have the same perfectly 
matched layers (PML) and the latter, as some post processes require this 
information to determine the integration order. 
If there exists a file 'pml.log' in the project directory, the PML will be 
based on the parameters in this file, otherwise, it will be created with the 
first scattering problem.
As the scattering problem has to be solved for many different frequencies,
the source file must contain the definition `Omega = %(omega)e`.

Custom parameters can be passed to the FEM solver via the property 'keys'
of Scattering. Please refer to the 
[parameter reference](https://www.docs.jcmwave.com/JCMsuite/html/ParameterReference/index.html)
and the documentation of the
[MATLAB interface](https://www.docs.jcmwave.com/JCMsuite/html/MatlabInterface/index.html)
of JCMsuite. Further details are provided in the examples.

### Adding quantities for expansion

For an efficient integration, not the total field is integrated along the
contours but the target quantities, which often are scalars. Arrays must 
be flattened to the size nx1 for integration and can be reshaped to their
original size later. Examples are the expansion of a cartesian export in 
'basicExample.m' and the expansion of the radiation pattern in 
'quadraticForms.m'. Both have user-defined plot functions. 
Please be aware that you may only expand holomorphic expressions.

Quantities are based on post processes that can be added with the method
'addPostProcess'. Several post processes are added automatically. Yet,
in order to make quantities available for expansion with the class 
'RieszProjection, the post processes must be adapted to the custom setup.
Therefore, quantities must be explicitly added to an instance of the class 
'Scattering' with the method 'addQuantity'. Modifications can be made 
by providing a custom keys structure for parameter substitution or by
changing the parameters used to compute the target from the returned 
values of the post process. The following examples, adapted from the script 
'metasurface.m', demonstrate how parameters and keys are used. 

```matlab
% Create interface to JCMsuite
sc = Scattering(projectFile, keys, wDir);
% The struct will be passed to the post process FourierTransform.jcmpt
ft_keys = struct('normalDirection','Z');
% The parameters defined here as key value pairs will be added as fields
% to the struct representing the quantity 'FourierTransform'
parameters = {'component',[1 0 0],'diffractionOrder',0};
sc.addQuantity('FourierTransform',ft_keys,parameters{:});
% now the quantity 'FourierTransform' is available for expansions
```

Above the quantity 'FourierTransform' has been adapted. In the next
example, the Fourier transform is used to define reflection and 
transmission. The interface to JCMsuite comes with the static method 'rt' 
which computes these quantities. In general, the results of post processes
must be evaluated with a function that takes the quantity (a struct with 
certain fields, e.g., 'jcmpt': path to the project file, 'resultbag':
a handle class containing results and 'parent': the instance of the class 
'Scattering') together with tags referring to results in the resultbag and
returns the quantity of interest. In this particular example the function 
'rt' expects the quantity to have a field 'nm' that provides a normalization
and is set with the method 'customizeQuantity'.

```matlab
% the power flux of the incoming plane wave for normalization
p = 0.5*keys.n_glass*sqrt(RieszProjection.eps0/RieszProjection.mu0);

keys_r = struct('normalDirection','Z'); % keys for reflection
keys_t = struct('normalDirection','-Z'); % keys for transmission

r_args = {keys_r,'evaluate',@Scattering.rt,'quadratic',true,'nm',p};
sc.customizeQuantity('Reflection','FourierTransform',r_args{:});
% The transmission is based on the reflection, i.e., everything is
% identical, except for the keys. 
sc.customizeQuantity('Transmission','Reflection',keys_t);
```

- Parameters that can be reset for any quantity: 
  - evaluate (function_handle) function used to evaluate the target from the 
    results of the post process
  - plot (function_handle) function used to plot the resulting expansion
    that will be called by RieszProjection/plot
  - quadratic (logical) quantities quadratic in E

The function handle evaluate takes the quantity, which is a struct, and
a cell array containing tags, each representing the result of a post
process applied to the solution of a scattering problem. It must return 
a double array of size (d,n) with d being a custom data dimension and 
n the number of tags. 

Currently, the following quantities are available for expansion: 

**Dipole emission and normalized decay rate of a point source**\
Usually, it is not necessary provide parameters. They are automatically
inferred from the keys structure for solving scattering problems which 
should have the corresponding fields.
  - Parameters:
    - [position](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/1b8a2f74fdcbfb512d7cda6d1e20efe5.html) (double) position of the point source
    - [strength](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/46e6d33084d3ba5c9f6e0f192b851c97.html) (double) strength and polarization of the source

**Electromagnetic field energy flux**\
The energy flux density is integrated across domain boundaries. 
  - Keys for parameter substitution:
    - [domainIdPairs](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/ac2c9ea9378bf25cb0457582e7c3c731.html) (double, optional) defines the interfaces
    - [interfaceType](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/47badf61609682dde1fbdcd3a81c84be.html) (char, default: 'ExteriorDomain') 

**Electromagnetic field absorption and electric field energy**\
Density integration provides the dissipated energy per second in absorptive
media.
  - Keys for parameter substitution: 
    - [domainIds](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/8f0ead8239d488ff0a815f372b46ca68.html) (double, default: 1) defines the domain of integration

**Farfield integral**\
If you add this quantity to the interface, it will add two quantities
for expansion to the class 'RieszProjection': 'Radiation' and 
'PhotonCollectionEfficiency'. The latter is available if for the
numerical aperture two values are provided, e.g., NA = [0.6 1], and
defined as the ratio of the energies emitted in the two solid angles
defined by the NAs. The NA is subdivided into intervals <= step and each 
interval is integrated with a 15-point Gauss-Kronrod quadrature rule.
  - Keys for parameter substitution: 
    - NA (double, default: [0.8 1]) numerical aperture
    - step (double, default: 10) step size in degrees
    - direction (char, default: 'up') choose direction ('up' or 'down')

**Radiation pattern**\
The radiation pattern displays the energy radiated into a given direction.
The evaluation points are defined by polar and azimuthal angle of spherical
coordinates. This quantity comes with a custom plot function: If two azimuthal 
angles are given a 2D representation is selected in polar coordinates and if
a full grid is defined a 3D visualisation is provided. 
 - Keys for parameter substiution: 
   - [gridPointsTheta](https://jcmwave.com/docs/ParameterReference/1f7ba10d1bf14a6d99efe55a9bc9b7d7.html) (double, default: 0:90) grid points of the polar angle
   - [gridPointsPhi](https://jcmwave.com/docs/ParameterReference/02f3e25777d49bf612ae117bcc5b547d.html) (double, default \[-90 90\]) grid points of the azimuthal angle

**Fourier transform**\
Fourier coefficients for discrete layouts that can, e.g., be used to
evaluate reflection and transmission coefficients. Alternatively, a single
component of a selected diffraction order can be expanded as a scalar 
quantity.
  - Keys for parameter substitution: 
    - [normalDirection](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/66663045236924c3d85488100e3b091a.html) (char, default: '-Y') defines the normal direction
  - Parameters
    - diffractionOrder (double) select diffraction order
    - component (double) select the component with a 1 by three vector

**Mode coupling**\
If the quantity of interest is e.g. the coupling to an optical fiber, the hidden
post process 'ModeCoupling.jcmpt' can be customized. A project directory
for solving the propagating mode problem has to be set up. Furthermore,
the keys structure passed to the method 'customizeQuantity' must contain a 
field with a corresponding '.jcmpt' file which contains the line: `Lambda0 
= %(lambda0)e`. As for the quadratic quantities independent solutions are
required, the corresponding polarizations must be aligned. For this, the 
electric field strength is exported. The position can be set with the field
'extract_polarization'. Furthermore, the position and orientation of the 
coordinate system used for the propagating modes relative to the one used 
for the scattering simulations must be specified. 
  - Keys for parameter substitution:
    - extract_polarization (double, default: [0 0 0]) position for export
    - [rotation](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/0796cb3ecafce08b7c1aceb3c493a3c2.html) (char, default: 'X:-Z:Y') orientation of the coord. system
    - [position](https://docs.jcmwave.com/JCMsuite/html/ParameterReference/86a5527157b0c62fb2777d0f35455d59.html?version=6.0.6) (double, default: [0 0 0]) position of the coord. system

**Cartesian Export**\
Sometimes it is of interest to get an impression of the field distribution in the
computational domain. For this purpose a cartesian export can be expanded and 
displayed with a custom plot function provided with the interface. 
  - Keys for parameter substitution
    - [GridPointsX](https://jcmwave.com/docs/ParameterReference/bd22d3c3a013c0aedf8aebf18b1ac8d9.html) (double, excludes NGridPointsX) grid points in x-direction
    - [NGridPointsX](https://jcmwave.com/docs/ParameterReference/63cfbf03388e1028794fa91ef0bd184c.html) (double, default: 100, excludes GridPointsX) number of equally distributed points in x-direction
    - [GridPointsY](https://jcmwave.com/docs/ParameterReference/293a54a8cc9f85802bcee2446d79270f.html) (double, excludes NGridPointsY) grid points in y-direction
    - [NGridPointsY](https://jcmwave.com/docs/ParameterReference/3b6725c3c0e940ba272111f6dcbda9e8.html) (double, default: 100, excludes GridPointsY) number of equally distributed points in y-direction
    - [GridPointsZ](https://jcmwave.com/docs/ParameterReference/b05053719407ef958e662759e565fa1d.html) (double, default: 0, excludes NGridPointsZ) grid points in z-direction
    - [NGridPointsZ](https://jcmwave.com/docs/ParameterReference/6f5ac313812b49c4e7460e9a9b315cb4.html) (double, excludes GridPointsZ) number of equally distributed points in z-direction
  - Parameters
    - component (char, default 'realx') defines the quantity to be displayed
The parameter 'component' can take the values real*, imag*, abs and log. The asterix 
stands for x, y and z. Be aware that the intensities will not sum up to the reference 
solution as taking the absolute value is not a holomorphic operation.


## Notes on a usage independent of JCMsuite

If you want to use your own software to solve the linear systems at the abscissae, 
you can still use the class 'RieszProjection'. The integration is based on 
quadrature rules of the form: 

```math
I(\omega_0) = \sum_i c_i(\omega_0) f(\omega_i)
```

Using the method 'getContours' you generate the weights $`c_i(\omega_0)`$
and the frequencies $`\omega_i`$
All you have to provide are the function values 
$`f(\omega_i)`$, e.g., solutions of the second-order Maxwell's equation:

```math
\nabla \times \mu^{-1}\nabla\times\bf{E}(\bf{r},\omega)-\omega^2\epsilon(\omega)\bf{E}(\bf{r},\omega)=i\omega\bf{J}(\bf{r}).
```

Usually, it is not necessary to integrate the total field. RPExpand is
designed to evaluate the target quantities at the integration points and, 
subsequently, to integrate them separately. As quantities based on quadratic 
forms require solutions of the two independent function values 
$`f(\omega_i)`$ and $`f(-\omega_i) = f^*(\omega_i^*)`$, the solution of the 
partial differential equation (the scattering problem) and the postprocess to 
extract the quantity of interest must be done in two independent steps. 

The interface to your favorite solver must be defined as a callable meeting
the following criteria:

The frequencies that discretize the contours are saved in a cell array 
'contours' with size (1, n) where n is the number of contours. Each element is a 
complex double array of shape (k,m), which corresponds to a contour. Here, k is the 
number of nodes per subinterval. Unless using higher-order quadrature methods, there 
is only one interval, i.e., m = 1 and k is the total number of integration points. 
In the first step, this cell array is passed to your function, which is expected to
return a cell array of size (1, n) or (n, 1) whose elements are objects of size (1, k\*m), 
e.g., a cell array containing numeric arrays, each of them representing the 
solution of a partial differential equation. In a second step, each element of 
the returned cell vector is passed to your custom function a second time with
an additional argument: `v = f(sc_results,quantity)`. The argument 'quantity' 
defines what has to be extracted from the results returned previously. It is 
defined as a 'char' vector (e.g., 'ElectricFieldEnergy'). The returned value 'v' 
is expected to be a cell containing numeric arrays of size (d, k\*m), where d is 
a custom number, which in many cases will be one, but can be larger as the examples 
for the expansion of the radiation pattern shown. In the case of a quantity quadratic 
in the solutions of the scattering problems, the first input is of size (n, 2) and the 
second column contains the solutions of the conjugated scattering problems. Otherwise, 
it is of shape (n, 1).
An example of a function meeting the described criteria is given in the file 
'quantumExample.m' which implements the quantum transition problem in 1D.

Eventually, you must provide a struct with the quantities available for expansion.
The field names must correspond to the names of the quantities and are structs
that must at least have the field 'quadratic' as it must be known beforehand if
solutions at complex conjugated frequencies are required. Furthermore, you can
add the field 'hidden' to specify if calling, e.g., the method 'expand' without 
arguments is supposed to list the corresponding quantity. 

## References

[[1]](http://dx.doi.org/10.1137/130931035) A. P. Austin et al.,
*SIAM J. Numer. Anal.* **52**, 1795 (2014).

[[2]](https://doi.org/10.1103/PhysRevA.98.043806) L. Zschiedrich et al.,
*Phys. Rev. A* **98**, 043806 (2018).

[[3]](http://dx.doi.org/%2010.1038/s42005-022-00977-1) F. Binkowski et al.,
*Commun. Phys.* **5**, 202 (2022).

[[4]](https://doi.org/10.48550/arXiv.2307.04654) F. Binkowski et al., 
arXiv:2307.04654 (physics.optics, 2023)

[[5]](https://doi.org/10.1103/PhysRevB.102.035432) F. Binkowski et al., 
*Phys. Rev. B* **102**, 035432 (2020)

[[6]](http://dx.doi.org/10.1364/JOSAA.428224) T. Wu et al.,
*J. Opt. Soc. Am. A* **38**, 1224 (2021).

[[7]](https://doi.org/10.1002/pssa.202200892) F. Betz et al.,
*Phys. Status Solidi A* **220**, 2200892 (2023)
