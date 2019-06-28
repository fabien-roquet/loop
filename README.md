# loop
Thermohaline loop, written in Fortran.


## Motivation

The thermohaline loop model consists of a fluid loop with infinitesimal section area along which a circulation is induced by applying fluxes of temperature and salt and/or by directly applying a torque.
The thermohaline model can serve as a conceptual model of the global ocean overturning circulation. It does not only find application in physical oceanography but
also in nuclear and solar energy engineering in the form of the so-called thermosyphons. 
Despite its simple setup, the thermohaline loop exhibits a variety of dynamical behaviors including instability and chaos. 
Apart from being interesting in its own right, the loop model has great educational value when it comes to the analysis of overturning flows. 

More information on the model equations can be found in the following publications:
* Pollmann, F., Roquet, F., and Madec, G., 2015. Effects of the Asymmetry between Surface and Interior Flow on the Dynamics of a Thermohaline Loop. Journal of Physical Oceanography. doi: 10.1175/JPO-D-15-0022.1
* Roquet F., Lindqvist R., Pollmann F., Ferreira D., Madec G., 2017. Stability of the thermohaline circulation examined with a one-dimensional fluid loop. Tellus A. doi: 10.1080/16000870.2017.1380490


## Versions

Two versions of the fluidloop are available. The first **fluidloop_v1** was used in Roquet et al. (2017). It is advised to use 
the more recent version **fluidloop_v2** which has been modified in a number of way:
1. Outputs are available in netcdf format, with machine precision. This format can be easily read with scientific programming tools
such as python or Matlab.
2. The diffusion and forcing are computed by default using an implicit scheme, which is unconditionally stable. This means that 
**fluidloop_v2** can be run with much greater time steps than **fluidloop_v1**.
3. An option is now available in **fluidloop_v2** to set variable diffusivities based on a constant energy of mixing scheme. This enables
studying the impact of variable stratification on the loop flow.


## Installation

Instructions for installation and execution of either **fluidloop_v1** or **fluidloop_v2** are available the
README.md files of their respective directory.


## Configuration files (for **fluidloop_v2**)

  1. The namelist is the configuration file that allows to configure the loop model. The structure of the namelist must be as follows:
  
    ```
    &namrun
      name_exp       = 'TEST'           ! name of experiment (used to name outputs). file max characters: 50
      nl             = 360              ! number of points around loop
      niter_t1       = 10000            ! number of iterations at t=1
      niter_max      = 10000  		  	  ! maximum number of timesteps
      niter_write    = 1000             ! number of iterations between writing output of output_1D
      use_netcdf     = .TRUE.           ! out put in netcdf or text file
      init_from_file = .TRUE.           ! start from init_state (T) or from rest (F)
      init_state = 'init_state.txt'     ! initial state for temperature, salinity and velocity. max characters: 50
    /

    &namconfig
      Zf = 0.5               	  	      ! height of source and sink of buoyancy -1<Zf<1
      R = 0.1                           ! inverse Rayleigh number
      lambda = 0.                     	! cabbeling parameter
      mu = 0.                           ! thermobaric parameter
      foldtrue = .TRUE.               	! Is the loop folded or not
      constant_diffusion = .FALSE.       ! Is the diffusion constant or not
      implicit_diffusion = .TRUE.       ! solve diffusion and forcing implicit (.FALSE. explicit)
    /

    &namforcing
      tau   = 0.                		    ! applied torque
      ref_t = 5.             	  	      ! Reference toward which temperature is relaxed
      xi_t  = 1.                        ! Temperature relaxation parameter
      F_s   = 0.15              	  	  ! Strength of fixed salinity forcing
    /
    ```

  2. Initial state
  
  If `init_from_file` in namelist is `.FALSE.`, all state variables are set to 0.
  
  Otherwise, initial values are read in the file named after the entry `init_state` in the namelist. The format is:
    - First line: initial velocity w
    - Line 2 to nl+1: index j, temperature t(j), salinity s(j)

  3. Forcing input
  
  Forcing parameters are read in the namelist and they are kept constant throughout the simulation.
  * tau    : applied torque
  * ref\_t : Reference temperature for relaxation (relaxed toward +ref\_t at source, and -ref\_t at sink)
  * xi\_t  : Temperature relaxation parameter
  * F\_s   : Magnitude of fixed salinity flux (+F\_s applied at source and -F\_s applied at sink)


## Output files

Three output files are created in the working directory, whose prefix `$EXP` is set by the `name_exp` entry in the `namrun` section of the namelist:

  1. `$EXP_config.txt` outputs the loop configuration, including namelist parameters and loop geometry factors
  
  2. `$EXP_out_0D.txt` provides snapshots of basic 0D properties of the loop at a frequency set by the `niter_write_0D` entry in the `namrun` section of the namelist:
    * niter
    * time
    * velocity
    * mass
    
  3. `$EXP_out_1D.txt` provides snapshots of basic 1D properties of the loop at a frequency set by the `niter_write_1D` entry in the `namrun` section of the namelist
    * niter: number of iteration
    * j: loop index (varies from 1 to nl)
    * theta: temperature
    * salt: salinity
    * sigma: density
    
  Examples of outputs are provided and can be used to test that the program has been properly installed. Note that the program execution will fail if output files already exist.

3. Three output files are created in the working directory, whose prefix `$EXP` is set by the `name_exp` entry in the `namrun` section of the namelist:
  1. `$EXP_config.txt` outputs the loop configuration, including namelist parameters and loop geometry factors
  2. `$EXP_out_0D.txt` provides snapshots of basic 0D properties of the loop at a frequency set by the `niter_write_0D` entry in the `namrun` section of the namelist
  3. `$EXP_out_1D.txt` provides snapshots of basic 1D properties of the loop at a frequency set by the `niter_write_1D` entry in the `namrun` section of the namelist
  More details on the format of these output files is provided in the section **Output files** below.
  
  Examples of outputs are provided and can be used to test that the program has been properly installed. Note that the program execution will fail if output files already exist.


## Matlab tools

Matlab tools are provided to:
* compute equilibrium tracer distribution when fixed-flux or relaxation forcings are applied: `loop_tracer_fixed.m` and `loop_tracer_relax.m`
* compute equilibrium state of the loop model with a semi-analytical value: `loop_equilibrium.m`
* analyze the linear stability of steady-states: `loop_jacobian.m`
* generate an initial-state file for the loop model: `loop_init_state.m`
* read model outputs in Matlab: `loop_read_out.m`

An example of use of Matlab tools is provided in the script `loop_test.m`

    
## Python tools

Python script `loop.py` is available thanks to the generous work of Julien Emile-Geay.
Also, a Jupyter notebook `plot_outputs_fluidloop.ipynb` developed by Eléa Lastes is also provided 
as an example of python scripts to analyse outputs of the model.

    
## Contributors

Fabien Roquet (Department of Marine Sciences at the University of Gothenburg). 
More information on the author on his [personal webpage](http://fabien-roquet.com/)

Other contributors: 
* Eléa Lastes (University of Marseille)
* Friederike Pollmann (University of Hamburg)
* Rickard Lindqvist (previously at Stockholm University)


## License

This software is distributed under the terms of the GNU GENERAL PUBLIC LICENSE v3. The GNU General Public License is a free, copyleft license for software and other kinds of works.

