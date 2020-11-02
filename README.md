# channel

> [Required dependencies](#dependencies)<br/>
> [Download and compile](#compile)<br/>
> [Preparing input files](#input)<br/>
> [Running](#running)

An exceptionally simple tool for Direct Numerical Simulation (DNS) of the incompressible Navier-Stokes equations 
in cartesian geometry, adapted from the engine by [Luchini & Quadrio, J. Comp. Phys. (2006)](https://www.sciencedirect.com/science/article/pii/S0021999105002871?via%3Dihub) and designed under the "Keep It Simple, Stupid" principle.

<img align="left" src="https://github.com/davecats/channel/blob/master/couette.png">
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


<a name="compile">

## Download and Compile

Once you have all required dependencies listed above, acquire the source code by cloning the git repository:

    git clone https://github.com/davecats/channel.git

Then compiling the code is as easy as hitting

    make
    
You may want to edit the *Makefile* and change the MPI/Fortran compiler, the optimization flag and the location of the *FFTW* library and headers.


<a name="input">

## Input files

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
```
- *nx* and *ny* are the number of modes in the statistically homogeneous x and z directions respectively. The corrisponding number of wavenumber is _2nx+1_ and _2nz+1_ (counting zero and doubling the numbers as there are both positive and negative wavenumbers); however, the actual number of x-modes stored in memory is _nx+1_ thanks to the Fourier transform of a real velocity field being Hermitian.
- *ny* is the number of points in the wall-normal y direction; the actual number of points, including walls, will be _ny+1_. However, the number of y points stored in memory is _ny+3_ due to the presence of ghost cells.
- *ni* is the scaling Reynolds number at which the simulation is performed.
- *ymin*, *ymax* specify the y coordinates of the walls; *a* is a parameter determining how points are distributed in the domain.

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=y_i=y_{min}+(y_{max}-y_{min})\left(\frac{tanh(a\,i/ny-1))}{tanh(a)}+1\right)" />
</p>

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


<a name="running">

## Running

The main program _channel_ must be run with mpi, in the following fashion:
```bash
mpirun -np number_of_cores /path/to/channel
```
where *number_of_cores* is indeed the number of cores used for parallel execution; it must be chosen so that _npxz_ is a divisor of _nx+1_ and _nzd_. While _nzd_ is printed out at the beginning of execution, _npxz_ can be found from *number_of_cores* with the following formula:
```
number_of_cores = npxz*npy
```
where _npy_ is defined in the file `mpi_transpose.f90` and is thus hardcoded. In a nutshell, two types of parallelisation are possible: in the statistically homogeneous directions x, z, and in the wall normal direction y. *npxz* specifies in how many subdivisions are performed for the parallelisation in x,z, while *npy* indicates the number of subdivisions in the y direction. Thus, the product of the two is the total number of parallel processes.

> Hint: *nzd* is always a power of 2 multiplied by 3; no other prime factors appear.

#### Quick note on needed input files

The directory in which _channel_ is called usually contains the following files:
- `dns.in` (always)
- `Dati.cart.out` (contains initial conditions; if absent, initial conditions are generated)
- `Runtimedata` (optional)

If *time_from_restart* is false, any *Runtimedata* file present in the directory where channel is called is overwritten. Otherwise, new timesteps will be appended at the end of file if it exists; a new one will be create if it doesn't.

## Contacts

Dr. Davide Gatti  
davide.gatti [at] kit.edu  
msc.davide.gatti [at] gmail.com  

Karlsruhe Institute of Technology  
Institute of Fluid Dynamics  
Kaiserstraße 10  
76131 Karlsruhe  
