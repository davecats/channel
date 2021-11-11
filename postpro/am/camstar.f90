
!--------------------------------------------------!
!-----------------     LEGEND     -----------------!
!--------------------------------------------------!

! mean(iv, iy)      where iv        U   W   Uy  Wy  Uyy Wyy P
!                                   1   2   3   4   5   6   7

! irs ----> index used for Reynolds stress
!   uu  vv  ww  uv  vw  uw
!   1   2   3   4   5   6

! uiujprofiles(irs, iterm, iy, largesmall)      where irs defined before, while iterm:
!   var prod    psdiss      ttrsp   tcross  vdiff   pstrain     ptrsp   PHIttrsp    PHIvdiff    PHIptrsp
!   1   2       3           4       5       6       7           8       9           10          11

! *VVd*(zf(iz), xf(ix), term, inout)
! : WARNING :  iz, ix indices here is scrambled as required by FFT -> use zf(iz), xf(ix)
! - inout = 1:5
! - term = 1:6
! inout=1 contains "inputs" (meaning, everything that is necessary for Fourier antitransform in order to go
! back to real domain). Products are calculated in real domain, and then transformed back to Fourier.
! All the "outputs" (everything after the Fourier antitransform, so basically products both in real and
! Fourier domains) are stored in inout=2.
! As for the "inputs" (inout=1):
! - term=1:4 are i,j velocity for large and small fields (both for ttrsp and tcross);
!   please use placeholders i_ft, j_ft summed to ilasm=1:2
! - term=6 is component k of the whole velocity field - not decomposed (only for tcross);
!   use placeholder uk_ft
! - term=5 is not used
! As for the "outputs" (inout=2):
! - for ttrsp, term=1:2 are large/small ui*uj product; use index ilasm
! - for tcross, term=1:4 are ui*uk, uj*uk for large/small;
!   please use placeholders i_ft, j_ft summed to ilasm=1:2
!   notice that two terms appear in tcross; for b=2, a=1:4,
!   the one starting with ui goes under i_fft, the one starting with uj goes under j_fft



program camstar
use dnsdata
implicit none

! stuff for parsing arguments

character(len=32) :: cmd_in_buf ! needed to parse arguments
character(len=32) :: arg

integer(C_INT) :: nfmin,nfmax,dnf,nftot
integer :: nmin,nmax,dn
integer :: nmin_cm = 0, nmax_cm = 0, dn_cm = 0

logical :: custom_mean = .FALSE.

! counters

integer :: ii ! counter initially used for parsing arguments
integer :: ix, iy, iz, iv ! for spatial directions and components

! global stuff

integer, parameter :: file_vel = 883
character(len=40) :: istring, currfname

real(C_DOUBLE), allocatable :: mean(:,:)

! camstar stuff

real :: retau, Lc, Lt, a0, acut
integer :: icut

! MPI stuff

TYPE(MPI_Datatype) :: mean_write_type, mean_inmem_type
TYPE(MPI_File) :: fh
integer :: ierror


! Program begins here
!----------------------------------------------------------------------------------------------------------------------------

    ! Init MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

    ! read arguments from command line
    call parse_args()
    if (custom_mean) then
        ! apply settings
        nfmin = nmin_cm; nfmax = nmax_cm; dnf = dn_cm
    end if
    nftot = ((nfmax - nfmin) / dnf) + 1
    
    ! setup
    call read_dnsin(.TRUE.) ! passing .TRUE. cause dns.in is in parent folder wrt cwd
    call camstar_setup()
    ! set npy to total number of processes
    ! this effectively DEACTIVATES PARALLELISATION IN X/Z
    ! and makes life easier (without significant losses in perf.)
    npy = nproc
    ! notice that this overrides value from dns.in
    call init_MPI(nx+1,nz,ny,nxd+1,nzd)
    call init_memory(.FALSE.) ! false flag avoids allocation of RHS related stuff!
    call init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB,.TRUE.) ! LAST FLAG (.TRUE.) SPECIFIES THAT REAL FFT TRANSFORM HAS AN ODD LOGICAL SIZE

    ! allocate stuff
    call init_camstar_mpitypes()
    allocate(mean(1:2, ny0-2:nyN+2))

    !---------------------------------------------------!
    !----------------- COMPUTE AVERAGE -----------------!
    !---------------------------------------------------!

    if (has_terminal) print *, "Computing average..."

    mean = 0
    do ii = nfmin, nfmax, dnf ! loop over files
        write(istring,*) ii
        open(file="../Dati.cart."//TRIM(ADJUSTL(istring))//".out", unit=file_vel, status="old", access="stream", action="read")
            do iy = ny0-2, nyN+2
                ! cumulate mean data
                mean(1,iy) = mean(1,iy) + real(fieldmap(file_vel, iy, 0, 0, 1)) ! U
                mean(2,iy) = mean(2,iy) + real(fieldmap(file_vel, iy, 0, 0, 3)) ! W
            end do
        close(file_vel)
    end do
    mean = mean / nftot ! divide by number of files

    ! revert to desired indices for calculation of tke
    nfmin = nmin; nfmax = nmax; dnf = dn
    nftot = ((nfmax - nfmin) / dnf) + 1

    !------------------------------------------------!
    !--------------- COMPUTE CAMSTAR ----------------!
    !------------------------------------------------!

    

    !---------------------------------------------------!
    !-----------------      WRITE      -----------------!
    !---------------------------------------------------!

    if (has_terminal) print *, "Saving to disk..."

    ! write mean data
    currfname = "CAMstar_mean.bin"
    call MPI_File_open(MPI_COMM_WORLD, trim(currfname), IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, fh)
        CALL MPI_File_set_view(fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION, mean_write_type, 'native', MPI_INFO_NULL)
        CALL MPI_File_write_all(fh, mean, 1, mean_inmem_type, MPI_STATUS_IGNORE)
    call MPI_File_close(fh)

    ! be polite and say goodbye
    if (has_terminal) print *, "Goodbye man!"

    !---------------------------------------------------!
    !-----------------     FINALISE    -----------------!
    !---------------------------------------------------!

    ! realease memory
    CALL free_fft()
    CALL free_memory(.FALSE.) 
    CALL MPI_Finalize()



contains !-------------------------------------------------------------------------------------------------------------------



    ! initialise largesmall stuff
    subroutine camstar_setup()

        ! read camstar.in
        open(15, file='camstar.in')
            read(15, *) retau, Lc
        close(15)

        ! calculate parameters
        Lt=retau*(2*PI/beta0); a0=2*PI/Lt; acut=2*PI/Lc; icut=FLOOR(acut/a0)

        ! write CAMstart.nfo
        if (has_terminal) then
            open(15, file='CAMstart.nfo')
                write(15,*) "This is camstar.f90."
                write(15,*) "nmin=", nmin, " nmax=", nmax, " dn=", dn, " nftot=", nftot
                write(15,*) "nmin_cm=", nmin_cm, " nmax_cm=", nmax_cm, " dn_cm=", dn_cm
                write(15,*) "Lt=", Lt, " Lc=", Lc, " icut=", icut, " 1000*2*PI/(beta0*icut)=", 1000*2*PI/(beta0*icut)
            close(15)
        end if

    end subroutine camstar_setup



    ! define and commit MPI filetypes
    subroutine init_camstar_mpitypes()

        ! define type for writing mean on disk
        CALL MPI_Type_create_subarray(2, [2, ny+3], [2, maxy-miny+1], [0,miny+1], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mean_write_type, ierror)
        CALL MPI_Type_commit(mean_write_type, ierror)
        CALL MPI_Type_create_subarray(2, [2, nyN-ny0+5], [2, maxy-miny+1], [0,miny-(ny0-2)], MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mean_inmem_type, ierror)
        CALL MPI_Type_commit(mean_inmem_type, ierror)

    end subroutine



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
            print *, "This is help for camstar.f90."
            print *, ""
            print *, "Calculates amplitude modulation map of channel flow."
            print *, "The simulation files will be looked for in the parent ('../)"
            print *, "of the working directory this program is called from."
            print *, "Usage:"
            print *, ""
            print *, "   mpirun -np no_proc camstar [-h] nfmin nfmax dn [-c --custom_mean] nmin_cm nmax_cm dn_cm"
            print *, ""
            print *, "A file 'camstar.in' must also be present in the working directory"
            print *, "when this program is called. See the template provided with"
            print *, "the source code."
            print *, ""
            print *, "Output will be written in the working directory."
            print *, "Output consists in CAMstar.bin and CAMstarxx.bin, where xx is"
            print *, "a number. The former contains averages, the latter istantaneous"
            print *, "values used to calculate the average. Both files are binary and"
            print *, "correspond to an array with indeces (iy2,iy1,i_cond,i_stat)."
            print *, "Indeces iy2 and iy1 are the position at which the large- and"
            print *, "small-scale signals are evaluated, respectively. Index i_cond"
            print *, "refers to conditioning: as per this program, 0 is a normal average,"
            print *, "1 is for positive and 2 for negative large-scale events."
            print *, "Index i_stat is instead referring to the statistic being calculated:"
            print *, "1 is the proper AM map, 2 and 3 are the autocorrelation of"
            print *, "small and large scales respectively."
            print *, ""
            print *, "Find more info on:"
            print *, "https://arxiv.org/abs/2109.09486"
            print *, ""
        end if
        call MPI_Finalize()
        stop
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
                case ('-c', '--custom_mean') 
                    custom_mean = .TRUE.
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmin_cm
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) nmax_cm
                    ii = ii + 1
                    call get_command_argument(ii, cmd_in_buf)
                    read(cmd_in_buf, *) dn_cm
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



end program