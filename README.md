# channel

> [Required dependencies](#dependencies)<br/>
> [Download and compile](#compile)<br/>
> [Preparing input files](#input)<br/>
> &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Important notice on restarting simulations](#notice_restart)<br/>
> [Parallelisation](#parallelisation)<br/>
> [Running](#running)<br/>
> [Runtime environment variables](#runtime_environment)<br/>
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
* *compact*: keeps the numerics in a small set of Fortran modules
* *parallel*: MPI x/z decomposition plus distributed wall-normal solves
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

```ini
[mesh]
nx = 191
ny = 384
nz = 189
alfa0 = 0.5
beta0 = 1.0
stretching = 1.5
ymin = 0.0
ymax = 2.0

[velocity]
ni = 12431.0
meanpx = 0.0
meanpz = 0.0
meanflowx = 0.0
meanflowz = 0.0
u0 = 0.0
uN = 0.0

[scalars]
nPhi = 0
meantx = 0.0
meantb = 0.0
t0 = 0.0
tN = 0.0

[timestepping]
deltat = 0.0
cflmax = 1.0
time = 0.0
dt_field = 30.0
dt_save = 30.0
t_max = 7000.0
time_from_restart = true
nstep = 999999
```
- *nx* and *nz* are the number of modes in the statistically homogeneous x and z directions respectively. The corrisponding number of points in physical space used for simulation is _2nx+1_ and _2nz+1_; however, this code is spectral, so x and z directions are in a spectral domain. Hence, the actual number of x-modes stored in memory is _nx+1_ thanks to the Fourier transform of a real velocity field being Hermitian. The number of z-modes is still _2nz+1_. See [domain](#domain).
- *ny* is the number of points in the wall-normal y direction; the actual number of points, including walls, will be _ny+1_. However, the number of y points stored in memory is _ny+3_ due to the presence of ghost cells. See [domain](#domain).
- *ni* is the scaling Reynolds number at which the simulation is performed.
- *ymin*, *ymax* specify the y coordinates of the walls; *a* is a parameter determining how points are distributed in the domain. See [domain](#domain).
- *meanpx* and *meanpz* prescribe a pressure gradient in the x and z directions; no pressure gradient is imposed if zero.
- *meanflowx* and *meanflowz* prescribe a flow rate in the x and z directions; no flow rate is imposed if zero.
- *u0* and *uN* represent a boundary condition; they are the x-component of the velocity at the walls.
- *time* prescribes the initial time instant for the simulation; such value is not used if *time_from_restart* is true.
- User can either specify a timestep _deltat_ or prescribe a maximum CFL (_cflmax_).
- *dt_field* specifies after how many time units a new snapshot is saved; the so saved snapshots can be used to calculate statistics.
- *dt_save* specifies after how many time units a restart file `Dati.cart.out` is generated. This __cannot__ be used to calculate statistics.
- *time_from_restart* is a boolean flag. If false, the restart file `Dati.cart.out` is used as the initial condition for the simulation, and the value *t0* is used as the initial value of time. If true, the initial value of time is read from the restart file.
- *t_max* and *nstep* specify respectively the final value of time and the maximum number of steps that one wants to achieve in a given run. After either of these two thresholds is reached, execution is terminated.

<a name="notice_restart">
 
### Important notice on restarting a simulation

When the program reaches the maximum time or the maximum number of steps (which are specified in dns.in), execution is interrupted and the last instant of time of the simulation is written to *Dati.cart.out*. So, if the previous simulation did reach maximum time or maximum iterations, the simulation can be restarted by just changing *time_from_restart* to .TRUE. in dns.in (also possibly increasing t_max in dns.in). The user does not have to manually modify *Dati.cart.out*.
 
However, if the previous run of the simulation did not reach either of the maximum time of the maximum iterations (for instance, if the simulation was stopped with CTRL+C or because of the wall-time limit on a cluster), the program will not update *Dati.cart.out* - which will remain the one of the previous simulation. 
So, if the user does not change *Dati.cart.out*, the new run will be in fact a repetition of the previous run. To avoid this, it is suggested that the user copies the last valid *Dati.cart.xx.out* to *Dati.cart.out* (making sure that such file is not corrupted, which is, that execution wasn't stopped while such file was being written). The simulation will start then from the time of *Dati.cart.xx.out*; we will call this time *Tx*.
 
Notice that the *Runtimedata* file will contain data about timesteps which come after *Tx* (because *Runtimedata* is written at each time step, and the previous run most likely performed quite some time steps after writing *Dati.cart.xx.out*). The program automatically recognises this: data at times greater than *Tx* are deleted from *Runtimedata*, and then the program simply continues to append to such file as usual.
  

 
<a name="parallelisation">

## Parallelisation

This program is parallelised with distributed memory. Each MPI rank belongs to one wall-normal group and one wall-parallel group. The number of subdivisions in the wall-parallel transform directions is stored in _npxz_; the number of subdivisions in the wall-normal direction is stored in _npy_. The total number of MPI ranks is
```
number_of_proc = npxz*npy
```
where _npy_ and _npxz_ are selected at runtime. By default, the MPI autotuner may choose _npxz_, _npy_ and the y-Schur hierarchy. Manual runs should set `CHANNEL_NPY`, `CHANNEL_NPXZ`, or both through the environment.

The x/z decomposition is used for FFT transposes. Because the transpose buffers use uniform counts, _npxz_ must divide both `nx+1` and `nzd = 3*nz`. The rank-local x/z block is
```
nx0 = ipxz*(nx+1)/npxz
nxN = (ipxz+1)*(nx+1)/npxz - 1
nz0 = ipxz*nzd/npxz
nzN = (ipxz+1)*nzd/npxz - 1
```
where `ipxz` is the rank coordinate inside the wall-parallel communicator.

The y decomposition splits the active unknown rows `1:ny-1` into contiguous slabs:
```
ny0 = 1 + ipy*(ny-1)/npy
nyN = (ipy+1)*(ny-1)/npy
```
where `ipy` is the rank coordinate inside the wall-normal communicator. Physical wall and ghost rows are only present on the end ranks, but the compact y solves also exchange the interface information needed by neighboring slabs.

### y-line Schur split

For each Fourier line, the compact wall-normal solve is a pentadiagonal linear system
```
A u = b.
```
With multiple y ranks, rank `r` owns a contiguous part of the line. Split its local unknowns into interior rows `i_r` and exposed interface rows `e_r`. The exposed rows are the two rows nearest each neighboring y slab; at physical walls the missing side is replaced by the wall boundary equations. In block form,
```
[ A_II  A_IE ] [ i_r ] = [ b_I ]
[ A_EI  A_EE ] [ e_r ]   [ b_E ].
```
The interior block is local to rank `r`, so it can be eliminated independently:
```
i_r = A_II^{-1}(b_I - A_IE e_r)
```
and therefore
```
(A_EE - A_EI A_II^{-1} A_IE) e_r = b_E - A_EI A_II^{-1} b_I.
```
This is the local Schur complement. The implementation packs this as up to four interface equations per rank and per Fourier line: one right-hand side plus coefficients for the four possible neighboring interface values. In other words each rank contributes an affine relation
```
e_r = c_r + G_r g_r,
```
where `g_r` are interface values owned by neighboring y slabs. Assembling all y ranks gives the reduced interface system
```
(I - G) e = c.
```
After this reduced system is solved, each rank reconstructs its eliminated interior rows using the first formula above.

The reduced system is solved hierarchically over the y communicator. A Schur pass list
```
p_1, p_2, ..., p_L
```
must satisfy
```
p_1*p_2*...*p_L = npy.
```
At level `ell`, groups of `p_ell` child systems are composed into a coarser Schur system. The span of a level is
```
P_ell = p_1*p_2*...*p_ell.
```
Ranks with the same `floor(ipy/P_ell)` are in the same parent group at that level. The final root system is solved, and interface values are propagated back down the same hierarchy. This keeps the global y solve distributed instead of gathering the full wall-normal line on one rank.

 
<a name="running">

## Running

The main program _channel_ must be run with mpi, in the following fashion:
```bash
mpirun -np number_of_proc /path/to/channel
```
where *number_of_proc* is indeed the number of processes used for parallel execution and must be specified by the user. If `CHANNEL_NPY` is set, the program calculates _npxz_ unless `CHANNEL_NPXZ` is also set. If neither variable is set, the autotuner can choose the decomposition; with autotuning disabled, the fallback is `npy = 1` and `npxz = number_of_proc`. The total number of processes must be chosen so that _npxz_ is a divisor of _nx+1_ and _nzd_; _nzd_ is printed out at the beginning of execution.

> Hint: *nzd* is always a power of 2 multiplied by 3; no other prime factors appear.

<a name="runtime_environment">

## Runtime environment variables

The following environment variables are optional. Boolean variables accept `1/0`, `true/false`, `yes/no`, and `on/off` forms unless stated otherwise.

| Variable | Default | Meaning |
| --- | --- | --- |
| `CHANNEL_NPY` | autotuned, otherwise `1` | Override the number of wall-normal MPI slabs. If `CHANNEL_NPXZ` is not set, `npxz = number_of_proc / CHANNEL_NPY`. |
| `CHANNEL_NPXZ` | autotuned, otherwise `number_of_proc` | Override the number of wall-parallel MPI groups. If `CHANNEL_NPY` is not set, `npy = number_of_proc / CHANNEL_NPXZ`. |
| `CHANNEL_MPI_AUTOTUNE` | enabled | `0`, `false`, `off`, or `no` disables autotuning. `report` prints a recommendation but does not apply it. When enabled and no manual decomposition variables are set, the autotuner may select _npxz_, _npy_, the Schur pass list, and the exchange mode. |
| `CHANNEL_MPI_AUTOTUNE_REPEATS` | `2` | Number of timed repeats per autotune candidate. The value is clamped to at least one. |
| `CHANNEL_Y_SCHUR_PASSES` | generated from _npy_ | Manual y-Schur pass list. Separators may be spaces, commas, semicolons, colons, `x`, or `X`. Each pass arity must be one of `2, 3, 4, 6, 8`, and the product must equal _npy_. |
| `CHANNEL_Y_SCHUR_EXCHANGE` | `auto` | Manual y-Schur exchange mode. Accepted values are `auto`, `default`, `alltoall`, `alltoallv`, `allgather`, and `allgatherv`. |
| `CHANNEL_Y_SCHUR_GLOBAL_EXCHANGE` | same as above | Backward-compatible alias checked before `CHANNEL_Y_SCHUR_EXCHANGE`. |
| `CHANNEL_OVERLAPPING` | disabled | Enables double-buffered overlap of x/z transpose communication and computation where supported. |
| `CHANNEL_DISABLE_RESTART_WRITE` | disabled | Skips writing the final `Dati.cart.out`, useful for benchmark/profiling runs that should not modify restart state. |
| `CHANNEL_YS_BATCH_MAX_COMPLEX` | `120000000` | Caps the complex workspace used by batched y-line endpoint solves. Lower it to reduce peak memory at the cost of smaller chunks. |
| `CHANNEL_YS_CHUNK_NX` | unlimited | Forces the maximum number of local x columns solved per y-line chunk. |
| `CHANNEL_YS_FORCE_CUSTOM_GPSV` | disabled | Nonzero integer forces the built-in pentadiagonal solver instead of the vendor sparse batched solver path. |

 
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

## Build instructions for Hunter

We use ftn here

- Build hipfort:  
  ```sh
  git clone git@github.com:ROCm/hipfort.git && cd hipfort
  mkdir build && cd build 
  cmake .. -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=~/hipfort-cray
  make -j # this takes a while (~30 mins)
  make install
  ```
- Build the channel project:
  ```
  mkdir build && cd build 
  cmake .. -DCMAKE_Fortran_COMPILER=ftn -Dhipfort_DIR=$HOME/hipfort-cray/lib/fortran/ftn/cmake/hipfort
  ```

Alternatively using AMD's new Fortran compiler:

- Download and extract the latest drop for RHEL on [the Radeon repo](https://repo.radeon.com/rocm/misc/flang/)
- You might need an older libffi for this to work. 
    - You can download the `x86_64` version from [Rocky Linux](https://dl.rockylinux.org/vault/centos/8.5.2111/BaseOS/x86_64/os/Packages/)
    - extract it using `rpm2cpio libffi-3.1-22.el8.x86_64.rpm | cpio -idmv` 
    - export LD_LIBRARY_PATH=$(pwd)/usr/lib64:$LD_LIBRARY_PATH
- Build using `export ROCM_AFAR_PATH=<...> && cmake .. -DCMAKE_Fortran_COMPILER=$ROCM_AFAR_PATH/bin/amdflang -DCMAKE_PREFIX_PATH=$ROCM_AFAR_PATH`


 
## Contacts

Dr. Davide Gatti  
davide.gatti [at] kit.edu  
msc.davide.gatti [at] gmail.com  

Karlsruhe Institute of Technology  
Institute of Fluid Dynamics  
Kaiserstraße 10  
76131 Karlsruhe  
