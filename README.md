# channel

> [Required dependencies](#dependencies)<br/>
> [Download and compile](#compile)<br/>
> [Preparing input files](#input)<br/>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Important notice on restarting simulations](#notice_restart)<br/>
> [Parallelisation](#parallelisation)<br/>
> [Running](#running)<br/>
> [Output files](#output)<br/>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Runtimedata](#notice_restart)<br/>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Velocity fields (Dati.cart*.out)](#velocity_fields)<br/>
> [Domain](#domain)<br/>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Why do we only store positive wavenumbers in the x direction?](#note_nxp1)<br/>
> [Postprocessing](#postpro)<br/>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Guidelines for usage/development](#standard_postprocessing)<br/>
> [Advanced: conditional compiler flags](#condcomp)<br/>

An exceptionally simple tool for Direct Numerical Simulation (DNS) of the incompressible Navier-Stokes equations 
in cartesian geometry, adapted from the engine by [Luchini & Quadrio, J. Comp. Phys. (2006)](https://www.sciencedirect.com/science/article/pii/S0021999105002871?via%3Dihub) and designed under the "Keep It Simple, Stupid" principle.

<img align="left" src="https://github.com/davecats/channel/blob/master/cover_pic.png">
<p> 
 <br/><br/><br/><br/><br/>
Turbulent Couette flow <br/> at a friction Reynolds number of Reτ=500 <br/> 3 752 787 600 DoF</p>
<br clear="left"/>
  
The code has been explicitly designed for shortness, compactness and simplicity, while still being parallel. Simplicity is preferred over excessive optimization. The code is optimized for reasonable parallel performance on up to O(2000) cores in 1024^3-sized problems. The main features are:

* *simple*: written with simplicity in mind 
* *compact*: consists of ~ 500 lines 
* *parallel*: two-dimensional pencil decomposition with MPI 
* *elegant*: data transposition with MPI interleaved datatypes and nonblocking communication
* *validated*: based on the engine developed by  [Luchini & Quadrio, J. Comp. Phys. (2006)](https://www.sciencedirect.com/science/article/pii/S0021999105002871?via%3Dihub)


<a name="dependencies">
 
## Required dependencies

* **GNU Make**
* **MPI**: version 3.1 or above with exposed mpi_f08 Fortran interface
* **FFTW**: version 3.x or above
* **FORTRAN**: any fortran f08 compliant compiler
* _Postprocessing only:_ the Python utilities included in this repository require the _channel_ python package ([andyandreolli/channel_pytools](https://github.com/andyandreolli/channel_pytools)). This can be installed with:
 ```bash
 pip install git+https://github.com/andyandreolli/channel_pytools#egg=channel
 ```


<a name="compile">

## Download and Compile

Once you have all required dependencies listed above, acquire the source code by cloning the git repository:

    git clone https://github.com/davecats/channel.git

Then, generate the `compiler.settings` file with the command:

    make configure

After adjusting this file as desired (e.g. with your compiler of choice and, where needed, the position of the FFTW library), compiling the code is as easy as hitting

    make
    
You can edit the *compiler.settings* also to use custom compiler flags. A preset of flags for debugging is also provided.


<a name="input">

## Input files

The directory in which _channel_ is called must contain the following files:

- `dns.in` (always)
- `Dati.cart.out` (contains initial conditions, always needed; if absent, initial conditions are generated)
- `Runtimedata` (optional; this is actually an output file)

We distinguish two use cases for this program. The user chooses one of these two cases with a boolean value *time_from_restart* in dns.in (see end of section).
1. A new simulation is started (*time_from_restart = .FALSE.*). In this case, the time of the simulation starts from zero. The initial condition is read from Dati.cart.out; if this file is absent, a initial condition is generated.
2. An old simulation is continued (*time_from_restart = .TRUE.*). The time of this new run starts from the last time of the previous run, which is read from Dati.cart.out. Also, the last velocity field of the previous run is read from Dati.cart.out, which will be used as the first velocity field of this new run.
 
In essence, **_Dati.cart.out_ always contains the initial field for the current run, also when an old simulation is continued**. As for *Runtimedata*, any *Runtimedata* file present in the directory where channel is called is overwritten, if *time_from_restart* is false (case 1). Otherwise (case 2, *time_from_restart = .TRUE.*), new timesteps will be appended at the end of the file if it exists; a new one will be create if it doesn't. Notice that Runtimedata is an output file containing simple statistics for each timestep; read more about it in the output section. **Please make sure to read the [notice about restarting a simulation](#notice_restart), as it contains useful information for both Dati.cart.out and Runtimedata.**
 
A file `dns.in` must be present in the directory `channel` is called from. Its structure needs to be something like this:

```FORTRAN
191 384 189                     ! nx, ny, nz
0.5d0 1.0d0                     ! alfa0 beta0  
12431.d0                        ! ni
1.5d0 0.0d0 2.0d0               ! a, ymin, ymax
.TRUE. 1 0.161436d0             ! CPI, CPItype, gamma
0.0d0  0.0d0                    ! meanpx, meanpz
0.0d0  0.0d0                    ! meanflowx meanflowz
0.0d0  0.0d0                    ! u0 uN
0.00d0 1.0d0 0.0d0              ! deltat, cflmax, t0
30.d0 30.d0 7000.d0 .TRUE.      ! dt_field, dt_save, t_max, time_from_restart
999999                          ! nstep
12                              ! npy
```
- *nx* and *nz* are the number of modes in the statistically homogeneous x and z directions respectively. The corrisponding number of points in physical space used for simulation is _2nx+1_ and _2nz+1_; however, this code is spectral, so x and z directions are in a spectral domain. Hence, the actual number of x-modes stored in memory is _nx+1_ thanks to the Fourier transform of a real velocity field being Hermitian. The number of z-modes is still _2nz+1_. See [domain](#domain).
- *ny* is the number of points in the wall-normal y direction; the actual number of points, including walls, will be _ny+1_. However, the number of y points stored in memory is _ny+3_ due to the presence of ghost cells. See [domain](#domain).
- *ni* is the scaling Reynolds number at which the simulation is performed.
- *ymin*, *ymax* specify the y coordinates of the walls; *a* is a parameter determining how points are distributed in the domain. See [domain](#domain).
- *CPI* is a flag that activates (if true) or deactivates constant-power-input-like forcing; this option has been designed for Poiseuille flows (with still walls), thus it will provide wrong results if the walls are moving (Couette flow). *CPItype* specifies the type of CPI-like forcing. _CPI = 0_ corresponds to a standard constant power input; in this case, *gamma* represents the fraction of power passed to the control, while the user should specify the desired power input by setting *ni* to be the power Reynolds number as in [here](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/global-energy-fluxes-in-turbulent-channels-with-flow-control/288CE28A721734161742427A0989E28D). Otherwise, _CPI = 1_ provides a constant ratio between laminar dissipation and total power input; the definition for laminar dissipation can still be found [here](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/global-energy-fluxes-in-turbulent-channels-with-flow-control/288CE28A721734161742427A0989E28D). In such case, the value of laminar dissipation is provided as variable *gamma*; the Reynolds number is arbitrary
- *meanpx* and *meanpz* prescribe a pressure gradient in the x and z directions; no pressure gradient is imposed if zero.
- *meanflowx* and *meanflowz* prescribe a flow rate in the x and z directions; no flow rate is imposed if zero.
- *u0* and *uN* represent a boundary condition; they are the x-component of the velocity at the walls.
- *t0* prescribes the initial time instant for the simulation; such value is not used if *time_from_restart* on next line is true.
- User can either specify a timestep _deltat_ or prescribe a maximum CFL (_cflmax_).
- *dt_field* specifies after how many time units a new snapshot is saved; the so saved snapshots can be used to calculate statistics.
- *dt_save* specifies after how many time units a restart file `Dati.cart.out` is generated. This __cannot__ be used to calculate statistics.
- *time_from_restart* is a boolean flag. If false, the restart file `Dati.cart.out` is used as the initial condition for the simulation, and the value *t0* is used as the initial value of time. If true, the initial value of time is read from the restart file.
- *tmax* and *nstep* specify respectively the final value of time and the maximum number of steps that one wants to achieve in a given run. After either of these two trhesholds is reached, execution is terminated.
- *npy* indicates in how many chunks the domain is divided in the y direction for parallelisation. See [parallelisation](#parallelisation).

<a name="notice_restart">
 
### Important notice on restarting a simulation

When the program reaches the maximum time or the maximum number of steps (which are specified in dns.in), execution is interrupted and the last instant of time of the simulation is written to *Dati.cart.out*. So, if the previous simulation did reach maximum time or maximum iterations, the simulation can be restarted by just changing *time_from_restart* to .TRUE. in dns.in (also possibly increasing t_max in dns.in). The user does not have to manually modify *Dati.cart.out*.
 
However, if the previous run of the simulation did not reach either of the maximum time of the maximum iterations (for instance, if the simulation was stopped with CTRL+C or because of the wall-time limit on a cluster), the program will not update *Dati.cart.out* - which will remain the one of the previous simulation. 
So, if the user does not change *Dati.cart.out*, the new run will be in fact a repetition of the previous run. To avoid this, it is suggested that the user copies the last valid *Dati.cart.xx.out* to *Dati.cart.out* (making sure that such file is not corrupted, which is, that execution wasn't stopped while such file was being written). The simulation will start then from the time of *Dati.cart.xx.out*; we will call this time *Tx*.
 
Notice that the *Runtimedata* file will contain data about timesteps which come after *Tx* (because *Runtimedata* is written at each time step, and the previous run most likely performed quite some time steps after writing *Dati.cart.xx.out*). The program automatically recognises this: data at times greater than *Tx* are deleted from *Runtimedata*, and then the program simply continues to append to such file as usual.
  

 
<a name="parallelisation">

## Parallelisation

This program is parallelised with distributed memory, meaning that computations are divided among different processes which communicate one with each other; each process has access to only a limited portion of data. More specifically, the simulation domain is divided into parts, each of which is given to a different process. Partitioning of the domain can be done:
1. in the statistically homogeneous x and z directions (both wall-parallel);
2. in the wall-normal direction y.
The number of subdivisions in the x/z directions is stored in variable _npxz_; the number of subdivisions in the wall-normal direction y is instead stored in _npy_. The total number of processes is thus:
```
number_of_proc = npxz*npy
```
where _npy_ is specified in _dns.in_, whereas npxz is automatically calculated from the number of processes. The number of processes is specified when the program is called. See [input files](#input) and [running](#running).

 
<a name="running">

## Running

The main program _channel_ must be run with mpi, in the following fashion:
```bash
mpirun -np number_of_proc /path/to/channel
```
where *number_of_proc* is indeed the number of processes used for parallel execution and must be specified by the user. The user thus specifies _npy_ and *number_of_proc*; the program thus calculates _npxz_ (see [parallelisation](#parallelisation) for more on *npxz* and *npy*). The total number of processes must be chosen so that _npxz_ is a divisor of _nx+1_ and _nzd_; _nzd_ is printed out at the beginning of execution.

> Hint: *nzd* is always a power of 2 multiplied by 3; no other prime factors appear.

 
<a name="output">
 
## Output files
 
The user starts the program channel in a generic directory; we will refer to this directory as _CWD_ (current working directory). This program stores all output files in CWD. These files are:
- the _Runtimedata_ file
- a series of _Dati.cart.ii.out_ files
- additional files, depending on conditional formatting (missing doc; for instance, immersed boundaries or body-forcing terms)
- a new version of _Dati.cart.out_ which overwrites the one used to start the simulation; keep in mind that this is actually an input file. The new _Dati.cart.out_ contains the last instant of time of the simulation, and gets written __only if maximum time (t_max) or maximum iterations (n_max) are reached__. Also read [this](#notice_restart).

<a name="runtimedata">

### Runtimedata

_Runtimedata_ is an ASCII file that gets written at each timestep; every line corresponds to an instant of time. Every column contains instead a different physical quantity. The column can be summarised as:
```
 time, dudy_bottom, dudy_top, dwdy_bottom, dwdy_top, fr_x, dpdx, fr_z, dpdz, XXX, deltat
```
Here we separated names by commas, but values in _Runtimedata_ are actually only separated by spaces/tabs.
- _time_ quite obviously is the simulation time; units inferred from dns.in.
- *dudy_xxx* and *dwdy_xxx* refer to the wall-normal gradients of the stream-wise (u) or span-wise (w) velocity components respectively; top and bottom correspond to the two different walls. Notice that data at the top wall is here changed in sign (so, at the top wall, -dudy is being written on Runtimedata).
- *fr_xxx* refers to the flow rate in the stream- (x) or span-wise direction.
- _dpdx_ and _dpdz_ are the pressure gradients in the stream- and span-wise directions respectively.
- XXX I don't know what this is, seriously. FIXME
- _deltat_ is the difference between the __next__ time (the time of the next line) and the current one.

<a name="velocity_fields">

### Velocity fields (Dati.cart*.out)

_Dati.cart.out_ and all the _Dati.cart.ii.out_ files are binary. They contain a header and a velocity field.
 
The header consists in 3 integers (_nx_, _ny_, _nz_) and seven double-precision floating point numbers (_alfa0_, _beta0_, _ni_, _a_, _ymin_, _ymax_, _time_). All of the variables dumped in the velocity-field-file are the same as dns.in, except for _time_, which indicates the simulation time of the velocity field being saved. Double-precision floating point numbers occupy 8 bytes each on disk; as for the integers, they usually take 4 bytes each. However, the size in bytes of an integer is not standardised, so it should be verified by the user for each machine/compiler. If indeed each integer takes 4 bytes on a given machine, the total length of the header is 68 bytes.
 
As for the velocity field, it is a 4-dimensional array of double-precision complex numbers. Each complex number thus occupies 16 bytes on disk. Please check section [domain](#domain) for info about the velocity field.
 

<a name="domain">
 
## Domain and variables
 
The simulation domain is three-dimensional; however, one dimension (the y-direction, wall-normal) refers to physical space, whereas the other two (x and z) refer to a Fourier domain. In other words, this code is spectral, as the unknowns of the simulation are Fourier-transformed in the stream- (x) and span-wise (z) directions. Consider for instance the velocity field; this is a four-dimensional array with indices:
```FORTRAN
V(iy,iz,ix,ic) ! FORTRAN order or column-major order 
```
The above is meant as FORTRAN ordering of the indices, meaning that _iy_ is the index that changes the fastest in memory. __BE CAREFUL__: if you are using C, or any other language that uses row-major order of the indices, the order of indices must be reversed, namely `V(ic,ix,iz,iy)`. The _iy_ index refers to the wall-normal (y) position in physical space, ix, iz refer to the stream- (x) and span-wise (z) Fourier modes.

1. Index _ix_ has dimension `nx+1` and bounds (0,nx); bounds are inclusive. Let _kx_ be the wavenumber ("Fourier variable") in the x direction; then, `kx=alfa0*ix`. For _alfa0_, see the [dns.in](#input).
2. Index _iz_ has dimension `2*nz + 1` and bounds (-nz,nz); bounds are inclusive. Let _kz_ be the wavenumber ("Fourier variable") in the z direction; then, `kz=beta0*iz`. For _beta0_, see the [dns.in](#input).
3. Index _iy_ has dimension ny+3 and bounds (-1,ny+1); bounds are inclusive. Points `iy=-1` and `iy=ny+1` correspond to ghost cells, whereas `iy=0` and `iy=ny` are the two walls of the channel, located at positions ymin and ymax respectively. In general, if y is the wall-normal spatial coordinate, it holds:
```FORTRAN
y(iy) = ymin + 0.5*(ymax-ymin)*(tanh(a*(2*iy/ny-1)) / tanh(a) + 1)
```
where ymin, ymax, a are defined in [dns.in](#input).

> Disclaimer:<br>
> Bounds of indeces are custom-defined in this program. This means, that indeces of arrays do not always start from 1, as it is normal in FORTRAN; sometimes, they are redefined, so that indeces start from 0 or some negative integer.<br>
> Re-defining index bounds is not always possible; for instance, if you are writing a Python script (or a C program) that reads output from this program, always remember that indices start from zero. Practically speaking, if you consider index iz, index `iz=-nz` in FORTRAN will be `iz=0` in C/Python; index `iz=nz` in FORTRAN will correspond to index `iz=2*nx` in C/Python.

 
<a name="note_nxp1">

### Why do we only store positive wavenumbers in the x direction?
 
TODO

---
<a name="postpro">
 
## Post-processing

<a name="standard_postprocessing">

### General guidelines for usage/development
 
The following rules/guidelines apply to all postprocessing tools (at least, the ones written in FORTRAN).
- When running the postprocessing tool, the current working directory (CWD) should be a subfolder of the directory containing the simulation. That is, the parent folder of CWD needs to contain the `dns.in`, the `Dati.cart.*.out` files and `Runtimedata`.
- The above mentioned sub-folder contains output and inputs that are specific to the postprocessing executable (eg., a settings file that is only used by the postprocessing executable).
- Instructions on how to run the executable can be accessed by running it with a flag `-h`.
- All inputs and arguments must be written to a `.nfo` file.
- At the end of execution, a string "EXECUTION COMPLETED ON ..." is appended to the `.nfo` file.
- The executable should check presence of input file, and abort if some is missing.
- The executable can be compiled from the Makefile in the root `channel` folder.

<a name="condcomp">
 
## Advanced: conditional compiler flags

TODO

 
## Contacts

Dr. Davide Gatti  
davide.gatti [at] kit.edu  
msc.davide.gatti [at] gmail.com  

Karlsruhe Institute of Technology  
Institute of Fluid Dynamics  
Kaiserstraße 10  
76131 Karlsruhe  
