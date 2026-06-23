#include "header.h"

module byte_workspace
  use, intrinsic :: iso_c_binding
  implicit none
  private

  integer(C_SIZE_T), parameter :: WORKSPACE_ALIGNMENT = 256_C_SIZE_T

  integer(C_INT8_T), allocatable, target, save :: workspace_raw(:)
  integer(C_SIZE_T), save :: workspace_base_index = 1_C_SIZE_T
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
  logical, save :: workspace_mapped = .false.
#endif
  integer(C_SIZE_T), save :: workspace_capacity = 0_C_SIZE_T
  integer(C_SIZE_T), save :: workspace_high_water = 0_C_SIZE_T
  logical, save :: workspace_owned = .false.
  character(len=64), save :: workspace_owner = ""

  public :: workspace_reserve, workspace_request, workspace_release, workspace_slice
  public :: workspace_capacity_bytes, workspace_high_water_bytes
  public :: workspace_finalize, workspace_align_offset, WORKSPACE_ALIGNMENT

contains

  subroutine workspace_reserve(nbytes)
    integer(C_SIZE_T), intent(in) :: nbytes
    integer(C_SIZE_T) :: alloc_bytes

    if (workspace_owned) then
      print *, "workspace_reserve: workspace already owned by ", trim(workspace_owner)
      error stop
    end if
    if (nbytes <= workspace_capacity) then
      workspace_high_water = max(workspace_high_water, nbytes)
      return
    end if

    call release_storage()
    alloc_bytes = nbytes + WORKSPACE_ALIGNMENT
    allocate (workspace_raw(alloc_bytes))
    call choose_base_index()
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    !$omp target enter data map(alloc: workspace_raw)
    workspace_mapped = .true.
#endif
    workspace_capacity = nbytes
    workspace_high_water = max(workspace_high_water, nbytes)
  end subroutine workspace_reserve

  subroutine workspace_request(nbytes, owner, ptr)
    integer(C_SIZE_T), intent(in) :: nbytes
    character(len=*), intent(in) :: owner
    type(C_PTR), intent(out) :: ptr

    if (workspace_owned) then
      print *, "workspace_request: workspace already owned by ", trim(workspace_owner), &
        "; requested by ", trim(owner)
      error stop
    end if

    call workspace_reserve(nbytes)
    workspace_high_water = max(workspace_high_water, nbytes)
    workspace_owned = .true.
    workspace_owner = owner
    ptr = c_loc(workspace_raw(workspace_base_index))
  end subroutine workspace_request

  subroutine workspace_release(owner)
    character(len=*), intent(in) :: owner

    if (.not. workspace_owned) then
      print *, "workspace_release: workspace is not owned; requested by ", trim(owner)
      error stop
    end if
    if (trim(owner) /= trim(workspace_owner)) then
      print *, "workspace_release: workspace owned by ", trim(workspace_owner), &
        "; release requested by ", trim(owner)
      error stop
    end if
    workspace_owned = .false.
    workspace_owner = ""
  end subroutine workspace_release

  subroutine workspace_slice(offset_bytes, ptr)
    integer(C_SIZE_T), intent(in) :: offset_bytes
    type(C_PTR), intent(out) :: ptr

    if (.not. workspace_owned) then
      print *, "workspace_slice: workspace is not owned"
      error stop
    end if
    if (offset_bytes > workspace_capacity) then
      print *, "workspace_slice: offset exceeds workspace capacity"
      error stop
    end if
    ptr = c_loc(workspace_raw(workspace_base_index + offset_bytes))
  end subroutine workspace_slice

  integer(C_SIZE_T) function workspace_capacity_bytes()
    workspace_capacity_bytes = workspace_capacity
  end function workspace_capacity_bytes

  integer(C_SIZE_T) function workspace_high_water_bytes()
    workspace_high_water_bytes = workspace_high_water
  end function workspace_high_water_bytes

  integer(C_SIZE_T) function workspace_align_offset(offset_bytes)
    integer(C_SIZE_T), intent(in) :: offset_bytes
    integer(C_SIZE_T) :: rem

    rem = modulo(offset_bytes, WORKSPACE_ALIGNMENT)
    if (rem == 0_C_SIZE_T) then
      workspace_align_offset = offset_bytes
    else
      workspace_align_offset = offset_bytes + WORKSPACE_ALIGNMENT - rem
    end if
  end function workspace_align_offset

  subroutine workspace_finalize()
    if (workspace_owned) then
      print *, "workspace_finalize: workspace still owned by ", trim(workspace_owner)
      error stop
    end if
    call release_storage()
    workspace_capacity = 0_C_SIZE_T
  end subroutine workspace_finalize

  subroutine choose_base_index()
    integer(C_INTPTR_T) :: addr
    integer(C_SIZE_T) :: rem

    addr = transfer(c_loc(workspace_raw(1)), addr)
    rem = modulo(int(addr, C_SIZE_T), WORKSPACE_ALIGNMENT)
    if (rem == 0_C_SIZE_T) then
      workspace_base_index = 1_C_SIZE_T
    else
      workspace_base_index = 1_C_SIZE_T + WORKSPACE_ALIGNMENT - rem
    end if
  end subroutine choose_base_index

  subroutine release_storage()
#if defined(HAVE_CUDA) || defined(HAVE_HIP)
    if (workspace_mapped) then
      !$omp target exit data map(delete: workspace_raw)
      workspace_mapped = .false.
    end if
#endif
    if (allocated(workspace_raw)) deallocate (workspace_raw)
    workspace_base_index = 1_C_SIZE_T
  end subroutine release_storage

end module byte_workspace
