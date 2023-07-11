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

#include "header.h"

MODULE mpi_transpose
  
  USE, intrinsic :: iso_c_binding
  USE, intrinsic :: iso_fortran_env
  USE typedef
  USE mpi_f08

  TYPE(MPI_Comm) :: MPI_CART_COMM, MPI_COMM_X, MPI_COMM_Y 
  integer(C_INT), save :: nproc,iproc,ierr,npx,ipx,npy,ipy,nzd,nx
  integer(C_INT), save :: nx0,nxN,nxB,nz0,nzN,nzB,block,nxpp
  integer(C_INT), save :: ny0,nyN,miny,maxy
  integer(C_INT), save :: TAG_LUDECOMP=100, TAG_LUDIVSTEP1=101, TAG_LUDIVSTEP2=102, TAG_DUDY=103
  complex(C_DOUBLE_COMPLEX), allocatable :: Ain(:),Aout(:)
  logical, save :: first,last,has_terminal,has_average
  TYPE(MPI_Datatype), save :: writeview_type,owned2write_type,vel_read_type,vel_field_type
#ifndef packunpack
  TYPE(MPI_Datatype), save :: Mdz,Mdx,cmpl
#endif
#ifdef nonblockingY
  ! Array of requests for nonblocking communication in linsolve
  TYPE(MPI_REQUEST), allocatable :: REQlinSolve(:), REQvetaTOuvw(:)
  ! Buffer arrays for manually buffered communication
  real(C_DOUBLE), allocatable       :: BUFlinSolve(:)
  real(C_DOUBLE), allocatable       :: BUFveta(:)
  integer(C_INT)                 :: szBUFlinSolve, elBUFlinSolve
  integer(C_INT)                 :: szBUFveta, elBUFveta
#endif

CONTAINS

! First we get the pack/unpack version of the transpose
! in blocking and nonblocking version
#ifdef packunpack
#ifndef nonblockingXZ
  !-------------- Transpose: Z to X --------------!
  !-----------------------------------------------!
  SUBROUTINE zTOx(Vz, Vx)  
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(1:,1:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:,1:)
    integer(C_SIZE_T) :: i,j
    integer(C_SIZE_T) :: ix,iz,one,num,zero,npxt
    integer(C_INT) :: ierr

    one = 1
    num = nzd/npx
    zero = 0
    npxt = npx

    Ain=0; Aout=0

    !Pack
    i=zero
    DO j=zero,npxt-one
      DO ix=nx0,nxN
        Ain(i:i+num-one)=Vz(one+j*num:(j+one)*num,ix-nx0+one)
        i=i+num
      END DO
    END DO
 
    !Send
    CALL MPI_Alltoall(Ain, nxB*nzB, MPI_DOUBLE_COMPLEX, Aout, nxB*nzB, MPI_DOUBLE_COMPLEX, MPI_COMM_X)

    !num = (nx+1)/npx
    !Unpack
    i=zero
    DO ix=0,nx
      Vx(ix+one,one:num)=Aout(i:i+num-one)
      i=i+num
    END DO
  END SUBROUTINE zTOx


  !-------------- Transpose: X to Z --------------!
  !-----------------------------------------------!
  SUBROUTINE xTOz(Vx,Vz)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:,1:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(1:,1:)
    integer(C_SIZE_T) :: i,j
    integer(C_SIZE_T) :: ix,iz,one,num,zero,npxt

    one = 1
    num = (nx+1)/npx
    zero = 0
    npxt = npx

    !Pack
    i=zero
    DO j=zero,npxt-one
      DO iz=nz0,nzN
        Ain(i:i+num-one)=Vx(one+j*num:(j+one)*num,iz-nz0+one)
        i=i+num
      END DO
    END DO

    !Send
    CALL MPI_Alltoall(Ain, nxB*nzB, MPI_DOUBLE_COMPLEX, Aout, nxB*nzB, MPI_DOUBLE_COMPLEX, MPI_COMM_X)

    !Unpack
    i=zero
    DO iz=0,nzd-1
      Vz(iz+one,one:num)=Aout(i:i+num-one)
      i=i+num
    END DO
  END SUBROUTINE xTOz
#else 
#error "No nonblockingXZ version of the pack/unpack transform yet!"
#endif

#else
! Then we get the datatype version of the transpose
! in blocking and nonblocking version
#ifndef nonblockingXZ
  !-------------- Transpose: Z to X --------------!
  !-----------------------------------------------!
  SUBROUTINE zTOx(Vz, Vx)  
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(1:,1:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:,1:)

    CALL MPI_Alltoall(Vz(:,:), 1, Mdz, Vx(:,:), 1, Mdx, MPI_COMM_X)

  END SUBROUTINE zTOx


  !-------------- Transpose: X to Z --------------!
  !-----------------------------------------------!
  SUBROUTINE xTOz(Vx,Vz)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:,1:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(1:,1:)

    CALL MPI_Alltoall(Vx(:,:), 1, Mdx, Vz(:,:), 1, Mdz, MPI_COMM_X)

  END SUBROUTINE xTOz
#else
  !-------------- Transpose: Z to X --------------!
  !-----------------------------------------------!
  SUBROUTINE zTOx(Vz, Vx)  
    complex(C_DOUBLE_COMPLEX), intent(in)  :: Vz(1:,1:)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vx(1:,1:)
    type(MPI_REQUEST), intent(inout) :: Rs

    CALL MPI_IAlltoall(Vz(:,:), 1, Mdz, Vx(:,:), 1, Mdx, MPI_COMM_X, Rs)

  END SUBROUTINE zTOx


  !-------------- Transpose: X to Z --------------!
  !-----------------------------------------------!
  SUBROUTINE xTOz(Vx,Vz)
    complex(C_DOUBLE_COMPLEX), intent(out) :: Vz(1:,1:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: Vx(1:,1:)
    type(MPI_REQUEST), intent(inout) :: Rs

    CALL MPI_IAlltoall(Vx(:,:), 1, Mdx, Vz(:,:), 1, Mdz, MPI_COMM_X, Rs)

  END SUBROUTINE xTOz
#endif
#endif


  !------- Divide the problem in 1D slices -------! 
  !-----------------------------------------------!
  SUBROUTINE init_MPI(nx,nz,ny,nxd,nzd)
    integer(C_INT), intent(in)  :: nx,nz,ny,nxd,nzd
    integer, parameter :: ndims = 4 
    integer :: array_of_sizes(ndims), array_of_subsizes(ndims), array_of_starts(ndims), ierror
    integer(C_INT) :: i,j,dims(1:2),coords(1:2)
    logical :: periods(1:2)=.false.,reorder=.true.
    TYPE(MPI_Datatype), save :: row, column, tmp
    integer(kind=MPI_ADDRESS_KIND) :: stride,lb
    integer(C_INT), allocatable, dimension(:) :: rankW, rankX, rankY
    nxpp=nx+1
    ! Define which process write on screen
    has_terminal=(iproc==0)
#ifdef nonblockingY
    TAG_LUDECOMP=0
    TAG_LUDIVSTEP1=(2*nz+1)
    TAG_LUDIVSTEP2=2*(2*nz+1)
    TAG_DUDY=3*(2*nz+1)
#endif
    ! Processor decomposition
#ifdef ycontiguous
    npx=nproc/npy;  
    dims(1)=npx; dims(2)=npy
    CALL MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, MPI_CART_COMM)
    CALL MPi_Cart_coords(MPI_CART_COMM, iproc, 2, coords)
    ipx=coords(1); ipy=coords(2)
    first=(ipy==0); last=(ipy==(npy-1))
    CALL MPI_Comm_split(MPI_CART_COMM, ipx, iproc, MPI_COMM_Y) ! Communicator with processes holding same X-slab
    CALL MPI_Comm_split(MPI_CART_COMM, ipy, iproc, MPI_COMM_X) ! Communicator with processes holding same Y-slab
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
    DO i=0,nproc-1
       IF (iproc==i) WRITE(*,*) "iproc=",iproc,"ipx=",ipx,"ipy=",ipy, "nx0=",nx0," nxN=",nxN," nxB=",nxB, "nz0=",nz0," nzN=",nzN," nzB=",nzB, "ny0=", ny0, "nyN=", nyN 
       CALL MPI_Barrier(MPI_COMM_WORLD)
    END DO
    FLUSH(output_unit)
#endif
#ifdef packunpack
    ! Allocate buffers for transposes
    ALLOCATE(Ain(0:block-1)); Ain=0
    ALLOCATE(Aout(0:block-1)); Aout=0
#else
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
#endif
    ! For READING VELOCITY, SETTING VIEW: datatype that maps velocity on disk to memory (it differs from writing: halo cells are read twice!)
    CALL MPI_Type_create_subarray(ndims, [ny+3, 2*nz+1, nxpp, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [ny0-1,0,nx0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_read_type, ierror)
    CALL MPI_Type_commit(vel_read_type, ierror)
    ! For READING VELOCITY, datatype in memory (avoids overflow) - NOTICE THAT THIS ARRAY IS FULL (NOT REALLY SUBARRAY)
    CALL MPI_Type_create_subarray(ndims, [nyN-ny0+5, 2*nz+1, nxB, 3], [nyN-ny0+5, 2*nz+1, nxB, 3], [0,0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, vel_field_type, ierror)
    CALL MPI_Type_commit(vel_field_type, ierror)
    ! For WRITING VELOCITY, SETTING VIEW: datatype to map distributed velocity array to file
    array_of_sizes = [ny+3, 2*nz+1, nxpp, 3] ! size along each dimension of the WHOLE array ON DISK
    array_of_subsizes = [maxy-miny+1, 2*nz+1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny+1, 0, nx0, 0] ! starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, writeview_type, ierror)
    CALL MPI_Type_commit(writeview_type, ierror)
    ! For WRITING VELOCITY, SKIPPING HALO CELLS: datatype with holes to skip halo cells and select only data to be written
    array_of_sizes = [(nyN+2)-(ny0-2)+1, 2*nz+1, nxB, 3] ! size along each dimension of the array IN MEMORY owned by each process
    array_of_subsizes = [maxy-miny+1, 2*nz+1, nxB, 3] ! size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    array_of_starts = [miny-(ny0-2), 0, 0, 0] ! starting position of each component; !!! IT'S ZERO BASED AND WRT TO ARRAY IN MEMORY !!!
    CALL MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes, array_of_starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_COMPLEX, owned2write_type, ierror)
    CALL MPI_Type_commit(owned2write_type, ierror)
#ifdef nonblockingY
    ! Allocate memory for send buffers [2*LU5Decomp (20 Doubles)  + 2*LeftLU5DivStep1 (4 Complex)]
    szBUFlinSolve = int(20*8+4*16+MPI_BSEND_OVERHEAD,C_SIZE_T)*int(2*nz+1,C_SIZE_T)
    elBUFlinSolve = szBUFlinSolve/8 + MERGE(1,0,MODULO(szBUFlinSolve,8)>0)
    ALLOCATE(BUFlinSolve(0:elBUFlinSolve-1))
    ! Allocate memory for send buffers [LeftLU5DivStep1 (2 Complex)]
    szBUFveta = int(2*16+MPI_BSEND_OVERHEAD,C_SIZE_T)*int(2*nz+1,C_SIZE_T)*int(nxN-nx0+1,C_SIZE_T)
    elBUFveta = szBUFveta/8 + MERGE(1,0,MODULO(szBUFveta,8)>0)
    ALLOCATE(BUFveta(0:elBUFveta-1))
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
