
!--------------------------------------------------!
!-----------------     LEGEND     -----------------!
!--------------------------------------------------!

! mean(iv, iy)      where iv        U   W
!                                   1   2

! irs ----> index used for Reynolds stress
!   uu  vv  ww  uv  vw  uw
!   1   2   3   4   5   6

#include "header.h"

program psd
use dnsdata
use ifport  ! this library is intel compiler specific; only used to create a directory (makedirqq)
            ! please comment this and substitute the call to makedirqq if using gfortran or other
implicit none

! stuff for parsing arguments

character(len=32) :: cmd_in_buf ! needed to parse arguments
character(len=32) :: arg

integer(C_INT) :: nfmin,nfmax,dnf,nftot
integer :: nmin,nmax,dn

integer :: wy_min, wy_max ! used to describe min and max iy of unique data of each process (includes wall, excludes ghost cells)
! notice: miny, maxy defined in mpi_transpose.f90 do the same, but they include ghost cells
! also: ny0, nyN defined in mpi_transpose.f90 do the same, but they exclude ALL halo cells (that means: both ghost and wall cells are excluded)

logical :: timestwo = .FALSE.

! counters

integer :: ii ! counter initially used for parsing arguments
integer :: ix, iy, iz ! for spatial directions and components
integer :: i, j, irs, c, iz_psd

! global stuff

integer, parameter :: file_vel = 883
character(len=40) :: istring, currfname

real(C_DOUBLE), allocatable :: mean(:,:), psddata(:,:,:,:)

character(40), parameter :: nfofile="psd.nfo" 

! MPI stuff

TYPE(MPI_Datatype) :: spectra_write_type, spectra_inmem_type
TYPE(MPI_File) :: fh
integer :: ierror
integer(MPI_OFFSET_KIND) :: offset

! Program begins here
!----------------------------------------------------------------------------------------------------------------------------

    ! Init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

    ! read arguments from command line
    call parse_args()
    call check_in_files()
    nftot = ((nfmax - nfmin) / dnf) + 1
    
    ! setup
    call read_dnsin(.TRUE.) ! passing .TRUE. cause dns.in is in parent folder wrt cwd
    ! set npy to total number of processes
    ! this effectively DEACTIVATES PARALLELISATION IN X/Z
    ! and makes life easier (without significant losses in perf.)
    npy = nproc
    ! notice that this overrides value from dns.in
    call init_MPI(nx+1,nz,ny,nxd+1,nzd)
    call init_memory(.FALSE.) ! false flag avoids allocation of RHS related stuff!
    ! final setup: define wy_min, wy_max
    wy_min = ny0; wy_max = nyN
    if (first) then; wy_min=wy_min-1; end if
    if (last) then; wy_max=wy_max+1; end if

    ! allocate stuff
    call init_mpitypes()
    allocate(mean(1:2,wy_min:wy_max))
    allocate(psddata(0:nx, 0:nz, wy_min:wy_max, 1:6))

    !---------------------------------------------------!
    !----------------- COMPUTE AVERAGE -----------------!
    !---------------------------------------------------!

    if (has_terminal) print *, "Computing average..."

    mean = 0
    do ii = nfmin, nfmax, dnf ! loop over files
        write(istring,*) ii
        open(file="../Dati.cart."//TRIM(ADJUSTL(istring))//".out", unit=file_vel, status="old", access="stream", action="read")
            do iy = wy_min, wy_max
                ! cumulate mean data
                mean(1,iy) = mean(1,iy) + real(fieldmap(file_vel, iy, 0, 0, 1)) ! U
                mean(2,iy) = mean(2,iy) + real(fieldmap(file_vel, iy, 0, 0, 3)) ! W
            end do
        close(file_vel)
    end do
    mean = mean / nftot ! divide by number of files

    !-------------------------------------------------!
    !--------------- COMPUTE SPECTRUM ----------------!
    !-------------------------------------------------!

    psddata = 0;

    do ii = nfmin, nfmax, dnf ! loop over files

        ! read velocity, pressure
        write(istring,*) ii
        currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
        call read_restart_file(currfname, V)

        ! remove average
        V(:,0,0,1) = V(:,0,0,1) - cmplx(mean(1,:))
        V(:,0,0,3) = V(:,0,0,3) - cmplx(mean(2,:))

#define u(cmp) V(iy,iz,ix,cmp)

        do iy = wy_min, wy_max
            do ix = nx0, nxN ! this is equivalent to nx = 0, nx owing to forcing parallelisation to be y-only
                c = 2; if (ix == 0) c = 1 ! multiplier for doubling points in x direction
                do iz = -nz, nz
                    do irs = 1, 6
                        call get_indexes(irs, i, j)
                        iz_psd = abs(iz)
                        psddata(ix, iz_psd, iy, irs) = psddata(ix, iz_psd, iy, irs) + c*cprod(u(i),u(j))
                    end do
                end do
            end do
        end do

    end do

    ! divide by nftot to obtain average
    psddata = psddata / nftot

    !---------------------------------------------------!
    !-----------------      WRITE      -----------------!
    !---------------------------------------------------!

    if (has_terminal) print *, "Saving to disk..."

    ! write metadata
    currfname = nfofile
    if (has_terminal) then
        open(15, file=currfname)
            write(15,*) "This is psd.f90."
            if (timestwo) write(15,*) "Correcting results by factor 2, to match older versions."
            write(15,*) ""
            write(15,*) "nmin", nmin 
            write(15,*) "nmax", nmax
            write(15,*) "dn", dn
            write(15,*) "nftot", nftot
            write(15,*) ""
        close(15)
    end if

    ! write to disk
    currfname = "psd.bin"
    call MPI_File_open(MPI_COMM_WORLD, trim(currfname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
        
        ! write psd data
        offset = 0
        CALL MPI_File_set_view(fh, offset, MPI_DOUBLE_PRECISION, spectra_write_type, 'native', MPI_INFO_NULL)
        CALL MPI_File_write_all(fh, psddata, 1, spectra_inmem_type, MPI_STATUS_IGNORE)

    call MPI_File_close(fh)

    if (has_terminal) call watermark()

    ! be polite and say goodbye
    if (has_terminal) print *, "Goodbye man!"

    !---------------------------------------------------!
    !-----------------     FINALISE    -----------------!
    !---------------------------------------------------!

    ! realease memory
    CALL free_memory(.FALSE.) 
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    ! define and commit MPI filetypes
    subroutine init_mpitypes()

        ! define type for writing uiujspectra on disk
        CALL MPI_Type_create_subarray(4, [nx+1, nz+1, ny+1, 6], [nx+1, nz+1, wy_max-wy_min+1, 6], [0,0,wy_min,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, spectra_write_type, ierror)
        CALL MPI_Type_commit(spectra_write_type, ierror)
        CALL MPI_Type_create_subarray(4, [nx+1, nz+1, wy_max-wy_min+1, 6], [nx+1, nz+1, wy_max-wy_min+1, 6], [0,0,0,0], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, spectra_inmem_type, ierror)
        CALL MPI_Type_commit(spectra_inmem_type, ierror)

    end subroutine



    subroutine get_indexes(irs, i, j)
    integer, intent(in) :: irs
    integer, intent (out) :: i, j
        select case(irs)
            case(1)
                i=1; j=1
            case(2)
                i=2; j=2
            case(3)
                i=3; j=3
            case(4)
                i=1; j=2
            case(5)
                i=2; j=3
            case(6)
                i=1; j=3
        end select
    end subroutine get_indexes



    function cprod(a, b) result(r)
    complex(C_DOUBLE_COMPLEX), intent(in) :: a, b
    real(C_DOUBLE) :: r
        r = dreal(dconjg(a)*b)
    end function cprod



!----------------------------------------------------------------------------------------------------------------------------
! less useful stuff here
!----------------------------------------------------------------------------------------------------------------------------



    subroutine print_help()
        if (iproc == 0) then
            print *, "Calculates (co-)spectra of all Reynolds stress tensor components in a channel flow."
            print *, "Statistics are calculated on files ranging from index nfmin to nfmax with step dn. Usage:"
            print *, ""
            print *, "   mpirun [mpi args] psd [-h -two] nfmin nfmax dn"
            print *, ""
            print *, "This program should be run so that dns.in is located in the parent folder."
            print *, "Results are output to psd.bin."
            print *, ""
            print *, "Flag '-two', or '--timestwo', is used to multiply results times two (for consistency with used framework)."
            call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
        end if
    end subroutine print_help



    ! read arguments
    subroutine parse_args()
        ii = 1
        do while (ii <= command_argument_count()) ! parse optional arguments
            call get_command_argument(ii, arg)
            select case (arg)
                case ('-h', '--help') ! call help
                    call print_help()
                    stop
                case ('-two', '--timestwo') ! multiply uiuj budget by two
                ! to match output of old version
                    timestwo = .TRUE.
                case default
                    if (command_argument_count() < 3) then ! handle exception: no input
                        print *, 'ERROR: please provide one input file as command line argument.'
                        stop
                    end if
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmin
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmax
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) dn
                    ! apply settings
                    nfmin = nmin; nfmax = nmax; dnf = dn
            end select
            ii = ii + 1
        end do
    end subroutine parse_args



    subroutine check_in_files()
    logical :: inputs_missing = .FALSE.

        if (iproc == 0) then

            print *
            print *, "Checking input files..."

            currfname = "../dns.in"
            call check_da_file(inputs_missing)

            ! check fields
            do ii = nfmin, nfmax, dnf ! loop over files
                ! check velocity
                write(istring,*) ii; currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
                call check_da_file(inputs_missing)
            end do

            ! return error if anything is missing
            if (inputs_missing) then
                print *
                print *, "ERROR: SOME INPUT FILES ARE MISSING"
                print *
                call MPI_Abort(MPI_COMM_WORLD, 0, ierr)
            end if

        end if

        call MPI_BARRIER(MPI_COMM_WORLD)

    end subroutine
    


    subroutine check_da_file(inputs_missing)
    logical, intent(inout) :: inputs_missing
    logical :: file_exists
    ! name of file to check needs to be written to currfname
        inquire(file=trim(currfname), EXIST=file_exists)
        if (.NOT. file_exists) then
            print *, "Missing ", trim(currfname)
            inputs_missing = .TRUE.
        end if
    end subroutine



    subroutine watermark()
    integer :: datetime(8)
    character(40) :: daystring, monthstring, yearstring, hstr, mstr, sstr

        call date_and_time(values=datetime)
        write(daystring,'(I0.2)') datetime(3)
        write(monthstring,'(I0.2)') datetime(2)
        write(yearstring,'(I0.4)') datetime(1)
        write(hstr,'(I0.2)') datetime(5)
        write(mstr,'(I0.2)') datetime(6)
        write(sstr,'(I0.2)') datetime(7)
        open(15, file=nfofile,access='append')
            write(15,*)
            write(15,*) "EXECUTION COMPLETED ON ", trim(daystring), "/", trim(monthstring), "/", trim(yearstring), " at ", trim(hstr), ":", trim(mstr), ":", trim(sstr) 
        close(15)

    end subroutine



end program