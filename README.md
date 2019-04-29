# kissDNS: channel-f90-mpi

An exceptionally simple tool for Direct Numerical Simulation (DNS) of the incompressible Navier-Stokes equations 
in cartesian geometry, adapted from the engine by [Luchini & Quadrio, J. Comp. Phys. (2006)](https://www.sciencedirect.com/science/article/pii/S0021999105002871?via%3Dihub) and designed under the "Keep It Simple, Stupid" principle.

The code has been explicitly designed for shortness, compactness and simplicity, while still being parallel. Simplicity is preferred over excessive optimization. The code is optimized for single-core performance and for reasonable parallel performance on up to O(100) cores in 512^3-sized problems. The main features are:

* *simple*: written with simplicity in mind 
* *compact*: consists of ~ 500 lines 
* *parallel*: with MPI 
* *elegant*: data transposition with MPI interleaved datatypes
* *validated*: based on the engine developed by  [Luchini & Quadrio, J. Comp. Phys. (2006)](https://www.sciencedirect.com/science/article/pii/S0021999105002871?via%3Dihub)
* *low-cost*: optimized for performing DNS on off-the-shelf hardware

