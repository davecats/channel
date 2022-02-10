program camstar
use dnsdata
implicit none

! stuff for parsing arguments

character(len=32) :: cmd_in_buf ! needed to parse arguments
character(len=32) :: arg

integer(C_INT) :: nfmin,nfmax,dnf,nftot
integer :: nmin,nmax,dn
integer :: nmin_cm = 0, nmax_cm = 0, dn_cm = 0

logical :: custom_mean = .FALSE., spatial_average = .FALSE.

! counters

integer :: ii ! counter initially used for parsing arguments
integer :: ix, iy, iz, iv, iy1, iy2 ! for spatial directions and components

! global stuff

integer, parameter :: file_vel = 883
character(len=40) :: istring, currfname
logical :: file_exists

integer :: nyh

real(C_DOUBLE), allocatable :: mean(:,:), CAM(:,:,:,:), CAMi(:,:,:,:), CAM0(:,:,:,:), mean_us2(:,:,:), mean_30(:,:,:), mean_ul(:,:,:), totmean_us2(:,:), totmean_us2_0(:,:)

! camstar stuff

real :: retau, Lc, Lt, a0, acut
integer :: icut

integer :: pm
integer, parameter :: d1=1, u1=2, d2=3, u2=4, large=1, small=2

character(40), parameter :: nfofile = 'CAMstart.nfo'

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
    call check_in_files()
    if (custom_mean) then
        ! apply settings
        nfmin = nmin_cm; nfmax = nmax_cm; dnf = dn_cm
    end if
    nftot = ((nfmax - nfmin) / dnf) + 1
    
    ! setup
    call read_dnsin(.TRUE.) ! passing .TRUE. cause dns.in is in parent folder wrt cwd
    call camstar_setup()
    ! set npy to 1
    ! this effectively DEACTIVATES PARALLELISATION IN Y
    ! and makes life easier (without significant losses in perf.)
    npy = 1
    ! notice that this overrides value from dns.in
    call init_MPI(nx+1,nz,ny,nxd+1,nzd)
    call init_memory(.FALSE.) ! false flag avoids allocation of RHS related stuff!
    call init_fft(VVdz,VVdx,rVVdx,nxd,nxB,nzd,nzB,.FALSE.,[4,2]) ! size of fft array is specified with optional flag

    ! allocate stuff
    call init_camstar_mpitypes()
    allocate(mean(1:2, ny0-2:nyN+2))
    nyh = ny/2 ! watch out: this is an integer division!
    allocate(CAM(1:3, 0:2, 0:nyh, 0:nyh))
    allocate(CAMi(1:3, 0:2, 0:nyh, 0:nyh))
    allocate(CAM0(1:3, 0:2, 0:nyh, 0:nyh))
    allocate(totmean_us2(0:nyh, 0:nyh))
    allocate(totmean_us2_0(0:nyh, 0:nyh))
    allocate(mean_us2(1:2, 0:nyh, 0:nyh))
    allocate(mean_ul(1:2, 0:nyh, 0:nyh))
    allocate(mean_30(1:2, 0:nyh, 0:nyh))

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
    
    if (has_terminal) print *, "Computing CAM..."

    ! set to zero arrays of interest
    CAM = 0 ! not necessary with CAMi: done for every field
    totmean_us2 = 0
    mean_us2 = 0
    mean_ul = 0

    ! remove old CAMstari.bin
    currfname = "CAMstari.bin"
    inquire(file=trim(currfname), EXIST=file_exists)
    if (file_exists) then
        open(unit=100, file=trim(currfname), status='old')
        close(100, status='delete')
    end if

    do ii = nfmin, nfmax, dnf ! loop over files

        ! read velocity
        write(istring,*) ii
        currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
        call read_restart_file(currfname, V)

        ! remove mean
        if (has_average) then
            if (spatial_average) then
                V(:,0,0,:) = 0
            else 
                V(:,0,0,1) = V(:,0,0,1) - cmplx(mean(1,:))
                V(:,0,0,3) = V(:,0,0,3) - cmplx(mean(2,:))
            end if
        end if

        CAMi = 0

        do iy2 = 0,nyh 
            do iy1 = 0,nyh

                ! do filtering at the 4 points of interest (y1 and y2 in both top and bottom halves of plane)
                call scale_filter(iy1, iy2)
                ! result of filtering is in rVVdx; see comments in subroutine scale_filter
                
                do iz = 1, nzB ! only own block in z
                    do ix = 1, (2*nxd) ! all points in x

                        CAMi(1,0,iy1,iy2) = CAMi(1,0,iy1,iy2) + (rVVdx(ix,iz,d1,small)**2)*rVVdx(ix,iz,d2,large) + (rVVdx(ix,iz,u1,small)**2)*rVVdx(ix,iz,u2,large)
                        CAMi(2,0,iy1,iy2) = CAMi(2,0,iy1,iy2) + (rVVdx(ix,iz,d1,small)**2)*(rVVdx(ix,iz,d2,small)**2) + (rVVdx(ix,iz,u1,small)**2)*(rVVdx(ix,iz,u2,small)**2)
                        CAMi(3,0,iy1,iy2) = CAMi(3,0,iy1,iy2) + rVVdx(ix,iz,d1,large)*rVVdx(ix,iz,d2,large) + rVVdx(ix,iz,u1,large)*rVVdx(ix,iz,u2,large)
                        totmean_us2(iy1,iy2) = totmean_us2(iy1,iy2) + (rVVdx(ix,iz,d1,small)**2) + (rVVdx(ix,iz,u1,small)**2)
                        
                        ! bottom wall, conditional
                        if (rVVdx(ix,iz,d2,large) > 0) then 
                            pm = 1
                        else
                            pm = 2
                        end if
                        CAMi(1,pm,iy1,iy2) = CAMi(1,pm,iy1,iy2) + (rVVdx(ix,iz,d1,small)**2)*rVVdx(ix,iz,d2,large)
                        CAMi(2,pm,iy1,iy2) = CAMi(2,pm,iy1,iy2) + (rVVdx(ix,iz,d1,small)**2)*(rVVdx(ix,iz,d2,small)**2)
                        CAMi(3,pm,iy1,iy2) = CAMi(3,pm,iy1,iy2) + rVVdx(ix,iz,d1,large)*rVVdx(ix,iz,d2,large)
                        mean_us2(pm,iy1,iy2) = mean_us2(pm,iy1,iy2) + (rVVdx(ix,iz,d1,small)**2)
                        mean_ul(pm,iy1,iy2) = mean_ul(pm,iy1,iy2) + rVVdx(ix,iz,d2,large)
                        
                        ! top wall, conditional
                        if (rVVdx(ix,iz,u2,large) > 0) then 
                            pm = 1
                        else
                            pm = 2
                        end if
                        CAMi(1,pm,iy1,iy2) = CAMi(1,pm,iy1,iy2) + (rVVdx(ix,iz,u1,small)**2)*rVVdx(ix,iz,u2,large)
                        CAMi(2,pm,iy1,iy2) = CAMi(2,pm,iy1,iy2) + (rVVdx(ix,iz,u1,small)**2)*(rVVdx(ix,iz,u2,small)**2)
                        CAMi(3,pm,iy1,iy2) = CAMi(3,pm,iy1,iy2) + rVVdx(ix,iz,u1,large)*rVVdx(ix,iz,u2,large)
                        mean_us2(pm,iy1,iy2) = mean_us2(pm,iy1,iy2) + (rVVdx(ix,iz,u1,small)**2)
                        mean_ul(pm,iy1,iy2) = mean_ul(pm,iy1,iy2) + rVVdx(ix,iz,u2,large)

                    end do
                end do                

            end do
        end do

        ! sum results of all processes into process 0 - CAMi only
        call MPI_REDUCE(CAMi, CAM0, 9*(nyh+1)*(nyh+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

        if (has_terminal) then

            ! divide by number of points in x and z (* 2 due to both half channels)
            CAM0 = CAM0 / (2*(2*nxd)*nzd)

            ! write camstari
            currfname = "CAMstari.bin"
            open(unit=100,file=trim(currfname),access="stream",action="write",position="append")
                write(100) CAM0
            close(100)

            ! cumulate CAM0 in CAM
            CAM = CAM + CAM0/nftot

        end if 

    end do

    ! sum results of all processes into process 0 - correction for conditional averages
    call MPI_REDUCE(mean_us2, mean_30, 2*(nyh+1)*(nyh+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    if (has_terminal) mean_us2 = mean_30
    call MPI_REDUCE(mean_ul, mean_30, 2*(nyh+1)*(nyh+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    if (has_terminal) mean_ul = mean_30
    call MPI_REDUCE(totmean_us2, totmean_us2_0, (nyh+1)*(nyh+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

    if (has_terminal) then
        ! apply correction for conditional averages
        mean_ul = mean_ul / ( (2*2*nxd*nzd) * nftot )
        mean_us2 = mean_us2 / ( (2*2*nxd*nzd) * nftot )
        totmean_us2_0 = totmean_us2_0 / ( (2*2*nxd*nzd) * nftot )
        do ii=1,2
            CAM(1,ii,:,:) = CAM(1,ii,:,:) - (mean_ul(ii,:,:) * totmean_us2_0(:,:))
        end do
    end if

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

    ! write CAMstar
    if (has_terminal) then
        currfname = "CAMstar.bin"
        open(unit=100,file=TRIM(currfname),access="stream",action="write")
            write(100) CAM
        close(100)
    end if

    ! write on nfo that execution has ended
    if (has_terminal) call watermark()

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
        if (iproc==0) then
            open(15, file=nfofile)
                write(15,*) "This is camstar.f90."
                write(15,*) ""
                if (spatial_average) then
                    write(15,*) "spatial_average = .TRUE."
                    write(15,*) "Removing spatial average to get fluctuation fields."
                    write(15,*) "That is, mode 0,0 of fluctuations is set to 0."
                else
                    write(15,*) "spatial_average = .FALSE."
                end if
                write(15,*) ""
                write(15,*) "nmin", nmin 
                write(15,*) "nmax", nmax
                write(15,*) "dn", dn
                write(15,*) "nftot", nftot
                write(15,*) ""
                write(15,*) "nmin_cm", nmin_cm
                write(15,*) "nmax_cm", nmax_cm
                write(15,*) "dn_cm", dn_cm
                write(15,*) "" 
                write(15,*) "Lt", Lt
                write(15,*) "Lc", Lc
                write(15,*) "icut", icut
                write(15,*) "1000*2*PI/(beta0*icut)", 1000*2*PI/(beta0*icut)
            close(15)
        end if

    end subroutine camstar_setup



    subroutine scale_filter(iy1, iy2)
    ! after this subroutine is called, rVVdx will contain the antitransforms
    ! of the u-field (1st component) in a xz plane at different positions, divided in
    ! large and small contributions. The structure of rVVdx is:
    ! rVVdx has size 2*(nxd+1),nzB,4,2
    ! 1. x-index: you actually need the first 2*nzd points, the last 3 are always 0
    ! 2. z-index
    ! 3. these points refer to different positions in y; use the parameters d1,u1,d2,u2,
    !    which correspond to y indeces y1, get_up(y1), y2, get_up(y2)
    ! 4. large or small field; see parameters large and small

        integer, intent(in) :: iy1, iy2
        integer :: up1, up2, icnt

        up1=get_up(iy1); up2=get_up(iy2) ! retrieves index of corresponding y location in upper half of channel

        ! small: icut+1:nz   ---   large: 0:icut

        ! copy V to VVdz array for FFT
        call copy4fft(d1,iy1,small); call copy4fft(u1,up1,small);  call copy4fft(d2,iy2,small); call copy4fft(u2,up2,small)
        call copy4fft(d1,iy1,large); call copy4fft(u1,up1,large);  call copy4fft(d2,iy2,large); call copy4fft(u2,up2,large)

        ! do fourier transform
        do icnt=1,2
            do iV=1,4
                call IFT(VVdz(:,:,iV,icnt)) ! first transform in z
                call MPI_Alltoall(VVdz(:,:,iv,icnt), 1, Mdz, VVdx(:,:,iv,icnt), 1, Mdx, MPI_COMM_X)
                VVdx(nx+2:nxd+1,1:nzB,iV,icnt)=0
                call RFT(VVdx(:,:,iV,icnt),rVVdx(:,:,iV,icnt)) ! second transform in x
            end do
        end do

    end subroutine



    subroutine copy4fft(iw,ir,lasm)
    integer, intent(in) :: iw, ir, lasm
    integer :: iwp1, iwp2, iwm1, iwm2, irp1, irp2, irm1, irm2

        ! WARNING: this differs from normal program, since there is no extension of the Fourier domain
        ! (or, no increase of spatial resolution to avoid dealiasing)
        ! i.e., no need to set anything zo zero in the arrays

        if (lasm == small) then

            ! iz >= 0: iz_fft = iz + 1
            iwp1 = (icut+2);                irp1 = (icut+1)
            iwp2 = (nz+1);                  irp2 = nz

            ! iz < 0: iz_fft = 2*(nz+1) + iz
            iwm1 = (nzd+1-nz);              irm1 = -nz
            iwm2 = (nzd+1-icut-1);          irm2 = (-icut-1)

        else

            ! iz >= 0: iz_fft = iz + 1
            iwp1 = 1;                       irp1 = 0
            iwp2 = (icut+1);                irp2 = icut

            ! iz < 0: iz_fft = 2*(nz+1) + iz
            iwm1 = (nzd+1-icut);            irm1 = (-icut)
            iwm2 = (nzd+1-1);               irm2 = (-1)

        end if

        ! actual copying
        VVdz(:,:,iw,lasm) = 0
        VVdz(iwp1:iwp2, :, iw, lasm) = V(ir, irp1:irp2, :, 1) ! iz >= 0
        VVdz(iwm1:iwm2, :, iw, lasm) = V(ir, irm1:irm2, :, 1) ! iz <= 0
    end subroutine



    function get_up(iy) result(iyup)
    integer, intent(in) :: iy
    integer :: iyup
        iyup = ny-iy
    end function



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
            print *, "   mpirun -np no_proc camstar [-h -s] nfmin nfmax dn [-c nmin_cm nmax_cm dn_cm]"
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
            print *, "OPTIONAL ARGUMENTS:"
            print *, ""
            print *, "   -s --spatial: sets mode 0,0 to 0 for fluctuations;"
            print *, "                 in other words, spatial average is removed"
            print *, "                 instead of spatial-temporal one."
            print *, "   -c --custom_mean: calculates mean on a different set of fields;"
            print *, "                     the fields on which mean is calculated are"
            print *, "                     indicated after the flag."
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
        integer :: ii
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
                case ('-s', '--spatial')
                    spatial_average = .TRUE.
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

            currfname = "camstar.in"
            call check_da_file(inputs_missing)

            ! check fields
            if (.NOT. custom_mean) then
                do ii = nfmin, nfmax, dnf ! loop over files
                    ! check velocity
                    write(istring,*) ii; currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
                    call check_da_file(inputs_missing)
                end do
            else ! if custom mean is used
                do ii = min(nmin,nmin_cm), max(nmax, nmax_cm), min(dn, dn_cm) ! loop over files
                    ! check velocity
                    write(istring,*) ii; currfname = trim("../Dati.cart."//TRIM(ADJUSTL(istring))//".out")
                    call check_da_file(inputs_missing)
                end do
            end if

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