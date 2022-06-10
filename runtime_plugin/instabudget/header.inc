    integer :: llz ! largest z index (iz) belonging to large scales -> large scales are -llz:llz
    real(C_DOUBLE) :: lc ! threshold wavelength for filter
    real(C_DOUBLE) :: ys, yl ! requested y positions for large and small scale signals
    logical :: remove_spatial_average ! flag to remove spatial or actual average from ul
    integer :: iys, iyl ! indexes of nearest y positions to the requested ones
    integer :: iv, ix, iz, is=1, il=1 ! a counter
    integer, parameter :: ft=777
    logical :: has_ys = .FALSE., has_yl = .FALSE.
    integer(MPI_OFFSET_KIND) :: no_reads = 0, fsize, disp
    real(C_DOUBLE) :: targtime

    type(MPI_COMM) :: MPI_COMM_L, MPI_COMM_S
    logical, allocatable :: gather_yl(:), gather_ys(:)
    integer :: snp=0, lnp=0 ! no of procs that have pieces of planes y=ys and y=yl
    integer, allocatable :: s_list(:), l_list(:) ! list of procs who have pieces of planes y=ys and y=yl
    type(MPI_GROUP) :: world_group, s_group, l_group
    
    complex(C_DOUBLE_COMPLEX), allocatable :: dVdy(:,:,:) ! array for gradient
#ifdef nonblockingY
    integer :: cnt_req
    type(C_PTR) :: ptrBUFveta
#endif

    type(MPI_File) :: fh_l, fh_us, fh_vs, fh_ws, fh_uvs, fh_duldy
    type(MPI_Datatype) :: yconst_plane_type
    integer(C_INT) :: ierror
    type(MPI_Status) :: status

#ifdef nonblockingXZ
     TYPE(MPI_REQUEST) :: req_rt(1:6)
#endif