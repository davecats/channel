! TODO: write into files as append + flag to remove average or 0 mode

#include "header.h"

module runtime
    use dnsdata
    use mpi_transpose
    use ffts
    implicit none
#include HEADER_RUNTIME
    
contains



    subroutine runtime_setup()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_SETUP_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_SETUP_SUBROUTINE   

    end subroutine



    subroutine runtime_finalise()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_FINALISE_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_FINALISE_SUBROUTINE

    end subroutine



    subroutine runtime_save()
    ! this subroutine is taken from a file in folder "runtime"
    ! the filename is contained in the macro RUNTIME_SAVE_SUBROUTINE
    ! such macro is set in header.h
#include RUNTIME_SAVE_SUBROUTINE

    end subroutine



#include RUNTIME_AUXILIARY_SUBROUTINES



end module