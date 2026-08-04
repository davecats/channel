#include "header.h"

! The computational grid and its distribution over ranks.
!
! These names used to be split between two modules for no reason other than
! history: the global sizes nx and nzd and every owned index range lived in
! mpi_transpose, while ny, nz, nxd and the spectral wavenumber arrays lived in
! dnsdata.  A reader had to know which half of the vocabulary lived where.
! They are all declared here now, and the modules that compute them fill them
! in: init_MPI (mpi_transpose) sets the decomposition, read_dnsin and
! init_memory (dnsdata) set the mesh and the wavenumbers.
!
! This module owns declarations only.  It deliberately depends on nothing, so
! it can sit underneath both the numerics and the physics.
module channel_grid

  use, intrinsic :: iso_c_binding

  implicit none

  !-------------------------------------------------------------------------
  ! Global mesh.
  !
  ! nx, ny, nz  : Fourier modes in x and z, intervals in y.
  ! nxd, nzd    : padded (3/2-rule) transform sizes.
  ! nPhi        : number of passive scalars.
  !-------------------------------------------------------------------------
  integer(C_INT), save :: nx, ny, nz
  integer(C_INT), save :: nxd, nzd
  integer(C_INT), save :: nPhi
  !$omp declare target(ny)

  !-------------------------------------------------------------------------
  ! Mesh definition: wall-parallel wavenumber spacing, wall-normal extent and
  ! the tanh stretching parameter.
  !-------------------------------------------------------------------------
  real(C_DOUBLE), save :: alfa0, beta0
  real(C_DOUBLE), save :: ymin, ymax, a

  !-------------------------------------------------------------------------
  ! Derived grid quantities.
  !
  ! y, dy    : wall-normal node positions and spacings.
  ! k2       : alfa^2 + beta^2 per (iz, ix) mode.
  ! ialfa    : i*alfa per ix,   ibeta : i*beta per iz.
  ! izd      : z mode -> padded transform index.
  ! dx, dz   : physical spacings used by the CFL estimate.
  ! factor   : 1/(2*nxd*nzd), the inverse-transform normalisation.
  !-------------------------------------------------------------------------
  real(C_DOUBLE), allocatable, save :: y(:), dy(:)
  real(C_DOUBLE), allocatable, save :: k2(:, :)
  complex(C_DOUBLE_COMPLEX), allocatable, save :: ialfa(:), ibeta(:)
  integer(C_INT), allocatable, save :: izd(:)
  real(C_DOUBLE), save :: dx, dz, factor

  !-------------------------------------------------------------------------
  ! Rank layout.  The ranks form an npxz-by-npy_grid grid; ipxz and ipy are
  ! this rank's position in it.
  !-------------------------------------------------------------------------
  integer(C_INT), save :: nproc, iproc
  integer(C_INT), save :: npy_grid = 1, npxz = 1, ipy = 0, ipxz = 0
  !$omp declare target(npy_grid, npxz, ipy, ipxz)

  !-------------------------------------------------------------------------
  ! Index ranges owned by this rank.  nxB and nzB are the corresponding block
  ! widths; the y range is ny0:nyN with two ghost rows on each side.
  !-------------------------------------------------------------------------
  integer(C_INT), save :: nx0, nxN, nxB
  integer(C_INT), save :: nz0, nzN, nzB
  integer(C_INT), save :: ny0, nyN
  !$omp declare target(nx0, nxN, nxB, nz0, nzN, nzB, ny0, nyN)

  !-------------------------------------------------------------------------
  ! y-Schur decomposition chosen at startup.  configure_mpi_decomposition sets
  ! all three (its arguments are intent(out)), so they need no initialiser --
  ! which also keeps this module free of any `use`.
  !-------------------------------------------------------------------------
  integer(C_INT), save :: npy = 1
  integer(C_INT), allocatable, save :: schur_pass_counts(:)
  integer(C_INT), save :: schur_exchange_mode

  !-------------------------------------------------------------------------
  ! has_terminal : this rank writes the human-readable diagnostics.
  ! has_average  : this rank owns the (0,0) mode, i.e. the mean profile.
  !-------------------------------------------------------------------------
  logical, save :: has_terminal, has_average

end module channel_grid
