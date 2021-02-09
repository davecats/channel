! Scalar SCAL was never assigned any value
! Element connectivity matrix is still empty
! one cool thing would be to parallelise the file format as well
! check out: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf

! also, one could switch to a structured grid

! also, a better argument parsing might be nice (like, taking mutltiple fields at once)


program out2vtk
! convert out file to VTK (XML) format for paraview
! syntax:
! mpirun -np 1 out2vtk Dati.cart.xx.out
use dnsdata
implicit none

character(len=32) :: cmd_in_buf ! needed to parse arguments

integer :: iv, ix, iz, iy, pos, ii, ierror
integer :: nxtot, nytot, nztot, ndim, nnos, nnoel, nel, nxloc, nyloc, nzloc, nnos_local, nel_local      
real*4, dimension(:,:), allocatable :: xyz, vec ! watch out for this type! must be the same as 
real*4, dimension(:), allocatable :: scal

integer :: zadd=1, yadd=1
integer, dimension(:), allocatable :: local_nnos_arr
integer :: start_node = 0, start_elem = 0

type(MPI_Datatype) :: wtype_3d, wtype_scalar, wtype_ien

integer, dimension(:,:), allocatable :: ien

integer :: etype

    ! read arguments
    if (COMMAND_ARGUMENT_COUNT() < 1) then
        print *, 'ERROR: please provide one input file as command line argument.'
        stop
    end if
    call get_command_argument(1, cmd_in_buf)
    read(cmd_in_buf, *) fname

    ! Init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    call read_dnsin()
    call init_MPI(nx+1,nz,ny,nxd+1,nzd)
    call init_memory(.TRUE.)
    CALL init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB)



    !------------!
    ! PARAMETERS !
    !------------!

    ndim  = 3  ! number of spatial dimension
    nnoel = 8  ! number of nodes per element
    etype = 11 ! element type: VTK_VOXEL (parallelepiped)

    ! GLOBAL QUANTITIES, which is, overall sizes of stuff
    nxtot = 2*nxd  ! no. points
    nytot = (ny + 1)   ! no. points
    nztot = nzd      ! no. points
    nnos  = nxtot*nytot*nztot             ! number of nodes
    nel   = (nxtot-1)*(nytot-1)*(nztot-1) ! number of elements

    ! things that are local to each process
    ! number of points in each direction
    nxloc = 2*nxd       ! these refer to UNIQUE data! It does not consider halo cells
    nyloc = min(ny,maxy)-max(0,miny)+1 ! these refer to UNIQUE data! It does not consider halo cells, nor ghost cells
    nzloc = nzB         ! these refer to UNIQUE data! It does not consider halo cells
    nnos_local = nxloc*nyloc*nzloc ! number of nodes OWNED BY PROCESS
    if (nz0 == 0) zadd = 0 ! this fixes number of elements to consider elements across chunks of data
    if (miny == -1) yadd = 0 ! this fixes number of elements to consider elements across chunks of data
    nel_local = (nxloc-1)*(nyloc-1+yadd)*(nzloc-1+zadd) ! number of elements OWNED BY PROCESS

    !-------------------!
    ! END OF PARAMETERS !
    !-------------------!



    ! allocate stuff
    allocate(xyz(ndim,nnos_local))    ! nodal coordinates
    allocate(vec(ndim,nnos_local))    ! vectorial field
    allocate(scal(nnos_local))        ! scalar field
    allocate(ien(nnoel,nel_local))    ! Element conectivity
    ien=0 ! FIXME: don't do this and actually calculate connectivity


    !---------------------------!
    ! MPI FILETYPES FOR WRITING !
    !---------------------------!

    ! each process tells the others how many local nodes it has
    allocate(local_nnos_arr(nproc))
    local_nnos_arr(iproc+1) = nnos_local
    do ii=0,(nproc-1)
        call MPI_BCAST(local_nnos_arr(ii+1),1,MPI_INTEGER,ii,MPI_COMM_WORLD)
    end do
    ! each process calculates its starting position in the global nodes list
    do ii=0,(iproc-1)
        start_node = start_node + local_nnos_arr(ii+1)
    end do
    ! same thing, but with elements
    local_nnos_arr = 0 ! recycle local_nnos_arr
    local_nnos_arr(iproc+1) = nel_local
    do ii=0,(nproc-1)
        call MPI_BCAST(local_nnos_arr(ii+1),1,MPI_INTEGER,ii,MPI_COMM_WORLD)
    end do
    ! each process calculates its starting position in the global nodes list
    do ii=0,(iproc-1)
        start_elem = start_elem + local_nnos_arr(ii+1)
    end do
    
    ! XYZ AND VEC - for writing (setting view)
    ! 0) number of dimensions of array you want to write
    ! 1) size along each dimension of the WHOLE array ON DISK
    ! 2) size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    ! 3) starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(2, [3, nnos], [3, nnos_local], [0, start_node], MPI_ORDER_FORTRAN, MPI_REAL, wtype_3d)
    !                             0)    1)         2)               3)
    CALL MPI_Type_commit(wtype_3d, ierror)

    ! SCALAR - for writing (setting view)
    ! 0) number of dimensions of array you want to write
    ! 1) size along each dimension of the WHOLE array ON DISK
    ! 2) size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    ! 3) starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(1, [nnos], [nnos_local], [start_node], MPI_ORDER_FORTRAN, MPI_REAL, wtype_scalar)
    !                             0)   1)         2)             3)
    CALL MPI_Type_commit(wtype_scalar, ierror)

    ! IEN - for writing (setting view)
    ! 0) number of dimensions of array you want to write
    ! 1) size along each dimension of the WHOLE array ON DISK
    ! 2) size of the PORTION of array TO BE WRITTEN BY EACH PROCESS
    ! 3) starting position of each component; !!! IT'S ZERO BASED !!!
    CALL MPI_Type_create_subarray(2, [nnoel, nel], [nnoel, nel_local], [0, start_elem], MPI_ORDER_FORTRAN, MPI_REAL, wtype_ien)
    !                             0)       1)             2)                 3)
    CALL MPI_Type_commit(wtype_ien, ierror)
    

    !----------------------------------!
    ! END OF MPI FILETYPES FOR WRITING !
    !----------------------------------!



    ! read file
    CALL read_restart_file(fname,V)

    ! for each y plane
    do iy = miny,maxy ! skips halo cells: only non-duplicated data
        
        ! skip ghost cells
        if (iy < 0) cycle
        if (iy > ny+1) cycle
        
        ! do fourier transform
        VVdz(1:nz+1,1:nxB,1:3,1)=V(iy,0:nz,nx0:nxN,1:3);         VVdz(nz+2:nzd-nz,1:nxB,1:3,1)=0;
        VVdz(nzd+1-nz:nzd,1:nxB,1:3,1)=V(iy,-nz:-1,nx0:nxN,1:3); 
        DO iV=1,3
            CALL IFT(VVdz(1:nzd,1:nxB,iV,1))  
            CALL MPI_Alltoall(VVdz(:,:,iV,1), 1, Mdz, VVdx(:,:,iV,1), 1, Mdx, MPI_COMM_X)
            VVdx(nx+2:nxd+1,1:nzB,iV,1)=0;    CALL RFT(VVdx(1:nxd+1,1:nzB,iV,1),rVVdx(1:2*nxd+2,1:nzB,iV,1))
        END DO
        ! rVVdx is now containing the antitransform
        ! it has size 2*(nxd+1),nzB,6,6
        ! but you only need (1:2*nxd,1:nzB,1:3,1)

        ! convert velocity vector for each process
        do ix=1,2*nxd
            do iz = 1,nzB
                pos = iz + (ix-1)*nzB + (iy-1)*nzB*2*nxd
                vec(:,pos) = rVVdx(ix,iz,:,1)
                xyz(1,pos) = (ix-1) * (2 * PI / alfa0)/(2*nxd-1)
                xyz(2,pos) = y(iy)
                xyz(3,pos) = (iz + nz0 - 2) * (2 * PI / beta0)/(nzd-1)
            end do
        end do

    end do
    
    call WriteXMLFormat(fname)

    ! realease memory
    CALL free_fft()
    CALL free_memory(.TRUE.) 
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    subroutine WriteXMLFormat(fname) ! xyz, ien, scal, vec, etype, ndim, nnoel, nnos, nel

        character(len = 40) :: fname

        character :: buffer*200, lf*1, offset*8, str1*8, str2*8
        integer   :: ivtk = 9, int, dotpos, i
        real*4    :: float ! this must have same type of xyz and vec and stuff
        integer(C_SIZE_T) :: nbytes_scal, nbytes_vec, nbytes_xyz, nbytes_ien, nbytes_offset, nbytes_etype
        integer(C_SIZE_T) :: ioff0, ioff1, ioff2, ioff3, ioff4, ioff5
        integer(C_SIZE_T) :: disp

        type(MPI_File) :: fh
        TYPE(MPI_Status) :: status

        lf = char(10) ! line feed character

        ! Layout for the Appended Data Section

        ! _ length | SCAL |      length | VEC |     lenght | XYZ |     lenght | IEN |     lenght | OFFSET |    lenght | ETYPE

        nbytes_scal   =         nnos * sizeof(float)
        nbytes_vec    = ndim  * nnos * sizeof(float)
        nbytes_xyz    = ndim  * nnos * sizeof(float)
        nbytes_ien    = nnoel * nel  * sizeof(  int)
        nbytes_offset =         nel  * sizeof(  int)
        nbytes_etype  =         nel  * sizeof(  int)

        ioff0 = 0                                   ! scal
        ioff1 = ioff0 + sizeof(int) + nbytes_scal   ! vec
        ioff2 = ioff1 + sizeof(int) + nbytes_vec    ! xyz
        ioff3 = ioff2 + sizeof(int) + nbytes_xyz    ! ien
        ioff4 = ioff3 + sizeof(int) + nbytes_ien    ! offset
        ioff5 = ioff4 + sizeof(int) + nbytes_offset ! etype

        ! adapt filename
        dotpos = scan(trim(fname),".", BACK= .true.)
        if ( dotpos > 0 ) fname = fname(1:dotpos)//"vtu"

        ! write xml
        if (has_terminal) then
            open(unit=ivtk,file=fname,form='binary')
                buffer = '<?xml version="1.0"?>'//lf                                                                                                  ; write(ivtk) trim(buffer)
                buffer = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf						  	                  ; write(ivtk) trim(buffer)
                buffer = '  <UnstructuredGrid>'//lf                                                                                                   ; write(ivtk) trim(buffer)
                write(str1(1:8),'(i8)') nnos
                write(str2(1:8),'(i8)') nel
                buffer = '    <Piece NumberOfPoints="'//str1//'" NumberOfCells="'//str2//'">'//lf                                                     ; write(ivtk) trim(buffer)
                buffer = '      <PointData> '//lf                                                                                                     ; write(ivtk) trim(buffer)
                write(offset(1:8),'(i8)') ioff0
                buffer = '         <DataArray type="Float32" Name="scalars" format="appended" offset="'//offset//'"           />'//lf                 ; write(ivtk) trim(buffer)
                write(offset(1:8),'(i8)') ioff1
                buffer = '         <DataArray type="Float32" Name="Velocity" NumberOfComponents="3" format="appended" offset="'//offset//'" />'//lf    ; write(ivtk) trim(buffer)
                buffer = '      </PointData>'//lf                                                                                                     ; write(ivtk) trim(buffer)
                buffer = '      <CellData>  </CellData>'//lf                                                                                          ; write(ivtk) trim(buffer)
                buffer = '      <Points>'//lf                                                                                                         ; write(ivtk) trim(buffer)
                write(offset(1:8),'(i8)') ioff2
                buffer = '        <DataArray type="Float32" Name="coordinates" NumberOfComponents="3" format="appended" offset="'//offset//'" />'//lf ; write(ivtk) trim(buffer)
                buffer = '      </Points>'//lf                                                                                                        ; write(ivtk) trim(buffer)
                buffer = '      <Cells>'//lf                                                                                                          ; write(ivtk) trim(buffer)
                write(offset(1:8),'(i8)') ioff3
                buffer = '        <DataArray type="Int32" Name="connectivity" format="appended" offset="'//offset//'" />'//lf                         ; write(ivtk) trim(buffer)
                write(offset(1:8),'(i8)') ioff4
                buffer = '        <DataArray type="Int32" Name="offsets" format="appended" offset="'//offset//'" />'//lf                              ; write(ivtk) trim(buffer)
                write(offset(1:8),'(i8)') ioff5
                buffer = '        <DataArray type="Int32" Name="types" format="appended" offset="'//offset//'" />'//lf                                ; write(ivtk) trim(buffer)
                buffer = '      </Cells>'//lf                                                                                                         ; write(ivtk) trim(buffer)
                buffer = '    </Piece>'//lf                                                                                                           ; write(ivtk) trim(buffer)
                buffer = '  </UnstructuredGrid>'//lf                                                                                                  ; write(ivtk) trim(buffer)
                buffer = '  <AppendedData encoding="raw">'//lf                                                                                        ; write(ivtk) trim(buffer)
                buffer = '_'                                                                                                                          ; write(ivtk) trim(buffer)
            close(ivtk)
        end if
        
        ! write scalar

        if (has_terminal) then
            open(unit=ivtk,file=fname,form='binary',position='append')
                write(ivtk) nbytes_scal
            close(ivtk)
        end if

        CALL MPI_Barrier(MPI_COMM_WORLD) ! barrier so every process retrieves the same filesize
        inquire(file=fname, size=disp) ! retrieve displacement

        CALL MPI_File_open(MPI_COMM_WORLD, TRIM(fname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
            CALL MPI_File_set_view(fh, disp, MPI_REAL, wtype_scalar, 'native', MPI_INFO_NULL)
            CALL MPI_File_write_all(fh, scal, size(scal), MPI_REAL, status)
        CALL MPI_File_close(fh)

        ! write vec

        if (has_terminal) then
            open(unit=ivtk,file=fname,form='binary',position='append')
                write(ivtk) nbytes_vec
            close(ivtk)
        end if
        
        CALL MPI_Barrier(MPI_COMM_WORLD) ! barrier so every process retrieves the same filesize
        inquire(file=fname, size=disp) ! retrieve displacement

        CALL MPI_File_open(MPI_COMM_WORLD, TRIM(fname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
            CALL MPI_File_set_view(fh, disp, MPI_REAL, wtype_3d, 'native', MPI_INFO_NULL)
            CALL MPI_File_write_all(fh, vec, size(vec), MPI_REAL, status)
        CALL MPI_File_close(fh)

        ! write xyz

        if (has_terminal) then
            open(unit=ivtk,file=fname,form='binary',position='append')
                write(ivtk) nbytes_xyz
            close(ivtk)
        end if

        CALL MPI_Barrier(MPI_COMM_WORLD) ! barrier so every process retrieves the same filesize
        inquire(file=fname, size=disp) ! retrieve displacement

        CALL MPI_File_open(MPI_COMM_WORLD, TRIM(fname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
            CALL MPI_File_set_view(fh, disp, MPI_REAL, wtype_3d, 'native', MPI_INFO_NULL)
            CALL MPI_File_write_all(fh, xyz, size(xyz), MPI_REAL, status)
        CALL MPI_File_close(fh)

        ! write element connectivity

        if (has_terminal) then
            open(unit=ivtk,file=fname,form='binary',position='append')
                write(ivtk) nbytes_ien
            close(ivtk)
        end if

        CALL MPI_Barrier(MPI_COMM_WORLD) ! barrier so every process retrieves the same filesize
        inquire(file=fname, size=disp) ! retrieve displacement

        CALL MPI_File_open(MPI_COMM_WORLD, TRIM(fname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
            CALL MPI_File_set_view(fh, disp, MPI_REAL, wtype_ien, 'native', MPI_INFO_NULL)
            CALL MPI_File_write_all(fh, ien, size(ien), MPI_REAL, status)
        CALL MPI_File_close(fh)

        if (has_terminal) then
            open(unit=ivtk,file=fname,form='binary',position='append')
                write(ivtk) nbytes_offset, (i,i=nnoel,nnoel*nel,nnoel)
                write(ivtk) nbytes_etype , (etype,i=1,nel)
                buffer = lf//'  </AppendedData>'//lf;           write(ivtk) trim(buffer)
                buffer = '</VTKFile>'//lf;                      write(ivtk) trim(buffer)
            close(ivtk)
        end if

    end subroutine



end program






