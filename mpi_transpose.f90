!============================================!
!                                            !
!            Distributed Traspose            !
!                  for the                   !
!      Direct Numerical Simulation (DNS)     !
!        of a turbulent channel flow         !
!                                            !
!============================================!
! 
! Author: Dr. Davide Gatti
! Date  : 15/Apr/2019
! 

#include 'header.h'

MODULE mpi_transpose
  
  USE, intrinsic :: iso_c_binding
  USE typedef
  USE mpi_f08

  TYPE(MPI_Comm) :: MPI_COMM_X, MPI_COMM_Y, MPI_COMM_NEXT, MPI_COMM_PREV
  integer(C_INT),save :: nproc,iproc,ierr,npx,ipx,npy=1,ipy
  integer(C_INT), save :: nx0,nxN,nxB,nz0,nzN,nzB,block
  integer(C_INT), save :: ny0,nyN,miny,maxy
  integer(C_INT), parameter :: TAG_LUDECOMP=100, TAG_LUDIVSTEP1=101, TAG_LUDIVSTEP2=102, TAG_DUDY=103
  logical, save :: first,last,has_terminal,has_average
  TYPE(MPI_Datatype), save :: Mdz,Mdx,cmpl,writeview_type,owned2write_type
#ifdef nonblockingY
  ! Array of requests for nonblocking communication in linsolve
  TYPE(MPI_REQUEST), allocatable :: REQlinSolve(:), REQvetaTOuvw(:)
#endif

CONTAINS

  !------- Divide the problem in 1D slices -------! 
  !-----------------------------------------------!
  SUBROUTINE init_MPI(nxpp,nz,ny,nxd,nzd)
    integer(C_INT), intent(in)  :: nxpp,nz,ny,nxd,nzd
    TYPE(MPI_Datatype) :: row,column,tmp
    integer(kind=MPI_ADDRESS_KIND) :: stride,lb
    integer, parameter :: ndims = 4 
    integer :: array_of_sizes(ndims), array_of_subsizes(ndims), array_of_starts(ndims), ierror
    ! Define which process write on screen
    has_terminal=(iproc==0)
    ! Processor decomposition
#ifdef ycontiguous
    npx=nproc/npy; ipx=iproc/npy; 
    CALL MPI_Comm_split(MPI_COMM_WORLD, ipx, iproc, MPI_COMM_Y) ! Communicator with processes holding same X-slab
    CALL MPI_Comm_rank(MPI_COMM_Y, ipy)
    first=(ipy==0); last=(ipy==(npy-1))
    CALL MPI_Comm_split(MPI_COMM_WORLD, ipy, iproc, MPI_COMM_X) ! Communicator with processes holding same Y-slab
#else
    npx=nproc/npy; ipy=iproc/npx; first=(ipy==0); last=(ipy==(npy-1));
    CALL MPI_Comm_split(MPI_COMM_WORLD, ipy, iproc, MPI_COMM_X) ! Communicator with processes holding same Y-slab
    CALL MPI_Comm_rank(MPI_COMM_X, ipx)
    CALL MPI_Comm_split(MPI_COMM_WORLD, ipx, iproc, MPI_COMM_Y) ! Communicator with processes holding same X-slab
#endif
    ! Calculate domain division in wall-normal direction
    ny0=1+ipy*(ny-1)/npy; nyN=(ipy+1)*(ny-1)/npy
    IF (first) THEN; miny=ny0-2; maxy=nyN;   ELSE; miny=ny0; END IF
    IF (last)  THEN; miny=ny0;   maxy=nyN+2; ELSE; maxy=nyN; END IF
    IF (first .AND. last) THEN;  miny=ny0-2; maxy=nyN+2;     END IF
    ! Calculate domain division 
    nx0=ipx*(nxpp)/npx;  nxN=(ipx+1)*(nxpp)/npx-1;  nxB=nxN-nx0+1;
    nz0=ipx*nzd/npx;   nzN=(ipx+1)*nzd/npx-1;   nzB=nzN-nz0+1;
    block=max(nxB*nzd,nxpp*nzB)
    has_average=(nx0==0)
    ! Communicators only with previous and next in Y
#ifdef mpiverbose
    WRITE(*,*) "iproc=",iproc,"ipx=",ipx,"ipy=",ipy, "nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN 
#endif
    ! MPI derived datatyped - basics
    CALL MPI_Type_contiguous(2,MPI_DOUBLE_PRECISION,cmpl)   !complex
    CALL MPI_Type_commit(cmpl)
    ! interlaved MPI datatypes - communicate plane of data
    CALL MPI_Type_vector(nxB,nzB,nzd,cmpl,row)
    lb=0; stride=8*2*nzB; CALL MPI_Type_create_resized(row,lb,stride,Mdz)
    CALL MPI_Type_commit(Mdz)    
    CALL MPI_Type_vector(nzB,1,nxd,cmpl,column)
    lb=0; stride=8*2;  CALL MPI_Type_create_resized(column,lb,stride,tmp)
    CALL MPI_Type_contiguous(nxB,tmp,Mdx)
    CALL MPI_Type_commit(Mdx)    
    ! For WRITING, SETTING VIEW: datatype to map distributed velocity array to file
    array_of_sizes = [ny+3, 2*nz+1, nxpp, 3] ! size along each dimension of the WHOLE array ON DISK
    array_of_subsizes = [maxy-miny+1, 2*nz+1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny+1, 0, nx0, 0] ! starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, writeview_type, ierror)
    CALL MPI_Type_commit(writeview_type, ierror)
    ! For WRITING, SKIPPING HALO CELLS: datatype with holes to skip halo cells and select only data to be written
    array_of_sizes = [(nyN+2)-(ny0-2)+1, 2*nz+1, nxB, 3] ! size along each dimension of the array IN MEMORY owned by each process
    array_of_subsizes = [maxy-miny+1, 2*nz+1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny-(ny0-2), 0, 0, 0] ! starting position of each component; !!! IT'S ZERO BASED AND WRT TO ARRAY IN MEMORY !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, owned2write_type, ierror)
    CALL MPI_Type_commit(owned2write_type, ierror)
#ifdef nonblockingY
    ! Allocate memory for requests
    ALLOCATE(REQlinSolve(1:(2*nz+1)*4), REQvetaTOuvw(1:(nxB+1)*(2*nz+1)))
#endif
  END SUBROUTINE init_MPI
  
#ifdef nonblockingY
  INCLUDE "rbparmat_nonblocking.f90"
#else
  INCLUDE "rbparmat_blocking.f90"
#endif

END MODULE mpi_transpose
