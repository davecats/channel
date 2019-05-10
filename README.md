# kissDNS: channel-f90-mpi

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

## Requisites 

* **MPI**: version 3.1 or above with exposed mpi_f08 Fortran interface
* **FFTW**: version 3.x or above
* **FORTRAN**: any fortran f08 compliant compiler
