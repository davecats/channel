#include "header.h"

program test_hipsparse_gpsv_allocations
  use, intrinsic :: iso_c_binding
#ifdef HAVE_HIP
  use hipfort_hipsparse_enums
#endif
  implicit none

#ifdef HAVE_HIP
  integer(C_INT), parameter :: M = 8_C_INT, BATCH_COUNT = 4_C_INT
  integer(C_INT), parameter :: HIP_MEMCPY_HOST_TO_DEVICE = 1_C_INT
  integer(C_INT), parameter :: HIP_MEMCPY_DEVICE_TO_HOST = 2_C_INT
  real(C_DOUBLE), parameter :: TOL = 1.0d-10

  complex(C_DOUBLE_COMPLEX), allocatable, target :: ds_host(:), dl_host(:), d_host(:), du_host(:), dw_host(:)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: x_rhs_host(:), x_exact_host(:), x_out_host(:)
  complex(C_DOUBLE_COMPLEX), allocatable, target :: ds_map(:), dl_map(:), d_map(:), du_map(:), dw_map(:), x_map(:)
  character(C_CHAR), allocatable, target :: buffer_map(:)
  type(C_PTR) :: handle
  type(C_PTR) :: ds_dev, dl_dev, d_dev, du_dev, dw_dev, x_dev, buffer_dev
  integer(C_SIZE_T) :: buffer_size
  integer(C_INT) :: status
  logical :: ok, all_ok

  interface
    function hipsparseCreate(handle) bind(c, name="hipsparseCreate")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipsparseCreate
      type(C_PTR) :: handle
    end function hipsparseCreate

    function hipsparseDestroy(handle) bind(c, name="hipsparseDestroy")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipsparseDestroy
      type(C_PTR), value :: handle
    end function hipsparseDestroy

    function hipsparseZgpsvInterleavedBatch_bufferSizeExt(handle, algo, m, ds, dl, d, du, dw, x, batch_count, pbuffer_size) &
      bind(c, name="hipsparseZgpsvInterleavedBatch_bufferSizeExt")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipsparseZgpsvInterleavedBatch_bufferSizeExt
      type(C_PTR), value :: handle
      integer(C_INT), value :: algo, m, batch_count
      type(C_PTR), value :: ds, dl, d, du, dw, x
      integer(C_SIZE_T) :: pbuffer_size
    end function hipsparseZgpsvInterleavedBatch_bufferSizeExt

    function hipsparseZgpsvInterleavedBatch(handle, algo, m, ds, dl, d, du, dw, x, batch_count, pbuffer) &
      bind(c, name="hipsparseZgpsvInterleavedBatch")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipsparseZgpsvInterleavedBatch
      type(C_PTR), value :: handle
      integer(C_INT), value :: algo, m, batch_count
      type(C_PTR), value :: ds, dl, d, du, dw, x, pbuffer
    end function hipsparseZgpsvInterleavedBatch

    function hipMalloc(ptr, size_bytes) bind(c, name="hipMalloc")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipMalloc
      type(C_PTR) :: ptr
      integer(C_SIZE_T), value :: size_bytes
    end function hipMalloc

    function hipFree(ptr) bind(c, name="hipFree")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipFree
      type(C_PTR), value :: ptr
    end function hipFree

    function hipMemcpy(dst, src, size_bytes, kind) bind(c, name="hipMemcpy")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipMemcpy
      type(C_PTR), value :: dst, src
      integer(C_SIZE_T), value :: size_bytes
      integer(C_INT), value :: kind
    end function hipMemcpy

    function hipDeviceSynchronize() bind(c, name="hipDeviceSynchronize")
      use, intrinsic :: iso_c_binding
      integer(C_INT) :: hipDeviceSynchronize
    end function hipDeviceSynchronize
  end interface

  allocate (ds_host(M*BATCH_COUNT), dl_host(M*BATCH_COUNT), d_host(M*BATCH_COUNT), du_host(M*BATCH_COUNT), dw_host(M*BATCH_COUNT))
  allocate (x_rhs_host(M*BATCH_COUNT), x_exact_host(M*BATCH_COUNT), x_out_host(M*BATCH_COUNT))
  call fill_pentadiagonal_system(ds_host, dl_host, d_host, du_host, dw_host, x_rhs_host, x_exact_host)

  status = hipsparseCreate(handle)
  call check_status(status, "hipsparseCreate")

  call query_buffer_size(handle, M, BATCH_COUNT, buffer_size)
  write (*, '(A,I0)') "hipSPARSE gpsv buffer size bytes: ", int(buffer_size, kind=kind(0))
  all_ok = .true.

  call run_hipmalloc_case(handle, ds_host, dl_host, d_host, du_host, dw_host, x_rhs_host, x_exact_host, buffer_size, ok)
  all_ok = all_ok .and. ok

  allocate (ds_map(M*BATCH_COUNT), dl_map(M*BATCH_COUNT), d_map(M*BATCH_COUNT), du_map(M*BATCH_COUNT), dw_map(M*BATCH_COUNT), x_map(M*BATCH_COUNT))
  ds_map = ds_host
  dl_map = dl_host
  d_map = d_host
  du_map = du_host
  dw_map = dw_host
  x_map = x_rhs_host

  allocate (buffer_map(max(1, int(buffer_size, kind=kind(1)))))

  !$omp target enter data map(to: ds_map, dl_map, d_map, du_map, dw_map, x_map) map(alloc: buffer_map)
  !$omp target data use_device_addr(ds_map, dl_map, d_map, du_map, dw_map, x_map, buffer_map)
  call run_case_from_ptrs("omp_map_arrays+omp_map_buffer", handle, c_loc(ds_map(1)), c_loc(dl_map(1)), c_loc(d_map(1)), &
                          c_loc(du_map(1)), c_loc(dw_map(1)), c_loc(x_map(1)), c_loc(buffer_map(1)), x_out_host, ok)
  !$omp end target data
  if (ok) ok = compare_solution("omp_map_arrays+omp_map_buffer", x_out_host, x_exact_host)
  all_ok = all_ok .and. ok

call copy_rhs_to_mapped_arrays(ds_map, dl_map, d_map, du_map, dw_map, x_map, ds_host, dl_host, d_host, du_host, dw_host, x_rhs_host)
  status = hipMalloc(buffer_dev, max(1_C_SIZE_T, buffer_size))
  call check_hip_status(status, "hipMalloc buffer_dev")
  !$omp target data use_device_addr(ds_map, dl_map, d_map, du_map, dw_map, x_map)
  call run_case_from_ptrs("omp_map_arrays+hipmalloc_buffer", handle, c_loc(ds_map(1)), c_loc(dl_map(1)), c_loc(d_map(1)), &
                          c_loc(du_map(1)), c_loc(dw_map(1)), c_loc(x_map(1)), buffer_dev, x_out_host, ok)
  !$omp end target data
  call check_hip_status(hipFree(buffer_dev), "hipFree buffer_dev")
  if (ok) ok = compare_solution("omp_map_arrays+hipmalloc_buffer", x_out_host, x_exact_host)
  all_ok = all_ok .and. ok

  status = hipMalloc(ds_dev, int(storage_size(ds_host(1))/8, kind=C_SIZE_T)*size(ds_host, kind=C_SIZE_T))
  call check_hip_status(status, "hipMalloc ds_dev")
  status = hipMalloc(dl_dev, int(storage_size(dl_host(1))/8, kind=C_SIZE_T)*size(dl_host, kind=C_SIZE_T))
  call check_hip_status(status, "hipMalloc dl_dev")
  status = hipMalloc(d_dev, int(storage_size(d_host(1))/8, kind=C_SIZE_T)*size(d_host, kind=C_SIZE_T))
  call check_hip_status(status, "hipMalloc d_dev")
  status = hipMalloc(du_dev, int(storage_size(du_host(1))/8, kind=C_SIZE_T)*size(du_host, kind=C_SIZE_T))
  call check_hip_status(status, "hipMalloc du_dev")
  status = hipMalloc(dw_dev, int(storage_size(dw_host(1))/8, kind=C_SIZE_T)*size(dw_host, kind=C_SIZE_T))
  call check_hip_status(status, "hipMalloc dw_dev")
  status = hipMalloc(x_dev, int(storage_size(x_rhs_host(1))/8, kind=C_SIZE_T)*size(x_rhs_host, kind=C_SIZE_T))
  call check_hip_status(status, "hipMalloc x_dev")

  call copy_host_to_device(ds_host, ds_dev)
  call copy_host_to_device(dl_host, dl_dev)
  call copy_host_to_device(d_host, d_dev)
  call copy_host_to_device(du_host, du_dev)
  call copy_host_to_device(dw_host, dw_dev)
  call copy_host_to_device(x_rhs_host, x_dev)

  !$omp target data use_device_addr(buffer_map)
  call run_case_from_ptrs("hipmalloc_arrays+omp_map_buffer", handle, ds_dev, dl_dev, d_dev, du_dev, dw_dev, x_dev, &
                          c_loc(buffer_map(1)), x_out_host, ok)
  !$omp end target data
  if (ok) ok = compare_solution("hipmalloc_arrays+omp_map_buffer", x_out_host, x_exact_host)
  all_ok = all_ok .and. ok

  call check_hip_status(hipFree(ds_dev), "hipFree ds_dev")
  call check_hip_status(hipFree(dl_dev), "hipFree dl_dev")
  call check_hip_status(hipFree(d_dev), "hipFree d_dev")
  call check_hip_status(hipFree(du_dev), "hipFree du_dev")
  call check_hip_status(hipFree(dw_dev), "hipFree dw_dev")
  call check_hip_status(hipFree(x_dev), "hipFree x_dev")

  !$omp target exit data map(delete: ds_map, dl_map, d_map, du_map, dw_map, x_map, buffer_map)
  deallocate (ds_map, dl_map, d_map, du_map, dw_map, x_map, buffer_map)

  status = hipsparseDestroy(handle)
  call check_status(status, "hipsparseDestroy")
  if (.not. all_ok) error stop "One or more hipSPARSE allocation cases failed"
  write (*, '(A)') "All hipSPARSE allocation cases passed."

#else
  write (*, '(A)') "Skipping: built without HAVE_HIP."
#endif

contains

#ifdef HAVE_HIP
  subroutine fill_pentadiagonal_system(ds, dl, d, du, dw, rhs, exact)
    complex(C_DOUBLE_COMPLEX), intent(out) :: ds(:), dl(:), d(:), du(:), dw(:), rhs(:), exact(:)
    integer(C_INT) :: i, b, p
    complex(C_DOUBLE_COMPLEX) :: sum

    ds = (0.0d0, 0.0d0)
    dl = (0.0d0, 0.0d0)
    d = (0.0d0, 0.0d0)
    du = (0.0d0, 0.0d0)
    dw = (0.0d0, 0.0d0)

    do i = 0, M - 1
      do b = 1, BATCH_COUNT
        p = i*BATCH_COUNT + b
        exact(p) = cmplx(1.0d0 + 0.1d0*i, -0.05d0*b, kind=C_DOUBLE)
        d(p) = cmplx(4.0d0 + 0.01d0*b, 0.0d0, kind=C_DOUBLE)
        if (i >= 1) dl(p) = cmplx(-0.75d0, 0.0d0, kind=C_DOUBLE)
        if (i >= 2) ds(p) = cmplx(-0.15d0, 0.0d0, kind=C_DOUBLE)
        if (i <= M - 2) du(p) = cmplx(-0.65d0, 0.0d0, kind=C_DOUBLE)
        if (i <= M - 3) dw(p) = cmplx(-0.10d0, 0.0d0, kind=C_DOUBLE)
      end do
    end do

    do i = 0, M - 1
      do b = 1, BATCH_COUNT
        p = i*BATCH_COUNT + b
        sum = d(p)*exact(p)
        if (i >= 1) sum = sum + dl(p)*exact((i - 1)*BATCH_COUNT + b)
        if (i >= 2) sum = sum + ds(p)*exact((i - 2)*BATCH_COUNT + b)
        if (i <= M - 2) sum = sum + du(p)*exact((i + 1)*BATCH_COUNT + b)
        if (i <= M - 3) sum = sum + dw(p)*exact((i + 2)*BATCH_COUNT + b)
        rhs(p) = sum
      end do
    end do
  end subroutine fill_pentadiagonal_system

  subroutine query_buffer_size(handle, m, batch_count, buffer_size)
    type(C_PTR), intent(in) :: handle
    integer(C_INT), intent(in) :: m, batch_count
    integer(C_SIZE_T), intent(out) :: buffer_size
    integer(C_INT) :: status
    complex(C_DOUBLE_COMPLEX), target :: dummy(1)

    dummy(1) = (0.0d0, 0.0d0)
    call check_hip_status(hipDeviceSynchronize(), "hipDeviceSynchronize before bufferSizeExt")
    status = hipsparseZgpsvInterleavedBatch_bufferSizeExt(handle, 0_C_INT, m, c_loc(dummy(1)), c_loc(dummy(1)), &
                                                          c_loc(dummy(1)), c_loc(dummy(1)), c_loc(dummy(1)), c_loc(dummy(1)), &
                                                          batch_count, buffer_size)
    call check_hip_status(hipDeviceSynchronize(), "hipDeviceSynchronize after bufferSizeExt")
    call check_status(status, "hipsparseZgpsvInterleavedBatch_bufferSizeExt")
  end subroutine query_buffer_size

  subroutine run_hipmalloc_case(handle, ds_h, dl_h, d_h, du_h, dw_h, x_h, exact_h, buffer_size, ok)
    type(C_PTR), intent(in) :: handle
    complex(C_DOUBLE_COMPLEX), intent(in), target :: ds_h(:), dl_h(:), d_h(:), du_h(:), dw_h(:), x_h(:), exact_h(:)
    integer(C_SIZE_T), intent(in) :: buffer_size
    logical, intent(out) :: ok
    type(C_PTR) :: ds_d, dl_d, d_d, du_d, dw_d, x_d, buffer_d
    integer(C_SIZE_T) :: bytes

    bytes = int(storage_size(ds_h(1))/8, kind=C_SIZE_T)*size(ds_h, kind=C_SIZE_T)
    call check_hip_status(hipMalloc(ds_d, bytes), "hipMalloc ds_d")
    call check_hip_status(hipMalloc(dl_d, bytes), "hipMalloc dl_d")
    call check_hip_status(hipMalloc(d_d, bytes), "hipMalloc d_d")
    call check_hip_status(hipMalloc(du_d, bytes), "hipMalloc du_d")
    call check_hip_status(hipMalloc(dw_d, bytes), "hipMalloc dw_d")
    call check_hip_status(hipMalloc(x_d, bytes), "hipMalloc x_d")
    call check_hip_status(hipMalloc(buffer_d, max(1_C_SIZE_T, buffer_size)), "hipMalloc buffer_d")

    call copy_host_to_device(ds_h, ds_d)
    call copy_host_to_device(dl_h, dl_d)
    call copy_host_to_device(d_h, d_d)
    call copy_host_to_device(du_h, du_d)
    call copy_host_to_device(dw_h, dw_d)
    call copy_host_to_device(x_h, x_d)

    call run_case_from_ptrs("hipmalloc_arrays+hipmalloc_buffer", handle, ds_d, dl_d, d_d, du_d, dw_d, x_d, buffer_d, x_out_host, ok)
    if (ok) ok = compare_solution("hipmalloc_arrays+hipmalloc_buffer", x_out_host, exact_h)

    call check_hip_status(hipFree(ds_d), "hipFree ds_d")
    call check_hip_status(hipFree(dl_d), "hipFree dl_d")
    call check_hip_status(hipFree(d_d), "hipFree d_d")
    call check_hip_status(hipFree(du_d), "hipFree du_d")
    call check_hip_status(hipFree(dw_d), "hipFree dw_d")
    call check_hip_status(hipFree(x_d), "hipFree x_d")
    call check_hip_status(hipFree(buffer_d), "hipFree buffer_d")
  end subroutine run_hipmalloc_case

  subroutine copy_rhs_to_mapped_arrays(ds_map, dl_map, d_map, du_map, dw_map, x_map, ds_h, dl_h, d_h, du_h, dw_h, x_h)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: ds_map(:), dl_map(:), d_map(:), du_map(:), dw_map(:), x_map(:)
    complex(C_DOUBLE_COMPLEX), intent(in) :: ds_h(:), dl_h(:), d_h(:), du_h(:), dw_h(:), x_h(:)
    ds_map = ds_h
    dl_map = dl_h
    d_map = d_h
    du_map = du_h
    dw_map = dw_h
    x_map = x_h
    !$omp target update to(ds_map, dl_map, d_map, du_map, dw_map, x_map)
  end subroutine copy_rhs_to_mapped_arrays

  subroutine run_case_from_ptrs(label, handle, ds_ptr, dl_ptr, d_ptr, du_ptr, dw_ptr, x_ptr, buffer_ptr, x_out, ok)
    character(*), intent(in) :: label
    type(C_PTR), intent(in) :: handle, ds_ptr, dl_ptr, d_ptr, du_ptr, dw_ptr, x_ptr, buffer_ptr
    complex(C_DOUBLE_COMPLEX), intent(out) :: x_out(:)
    logical, intent(out) :: ok
    integer(C_INT) :: status

    call check_hip_status(hipDeviceSynchronize(), "hipDeviceSynchronize before "//trim(label))
  status = hipsparseZgpsvInterleavedBatch(handle, 0_C_INT, M, ds_ptr, dl_ptr, d_ptr, du_ptr, dw_ptr, x_ptr, BATCH_COUNT, buffer_ptr)
    write (*, '(A,1X,A,1X,I0)') "solve status", trim(label)//":", status
    if (status /= HIPSPARSE_STATUS_SUCCESS) then
      ok = .false.
      return
    end if

    call check_hip_status(hipDeviceSynchronize(), "hipDeviceSynchronize after "//trim(label))
    call copy_device_to_host(x_ptr, x_out)
    ok = .true.
  end subroutine run_case_from_ptrs

  logical function compare_solution(label, actual, expected) result(matches)
    character(*), intent(in) :: label
    complex(C_DOUBLE_COMPLEX), intent(in) :: actual(:), expected(:)
    integer :: i
    real(C_DOUBLE) :: err

    matches = .true.
    do i = 1, size(actual)
      err = abs(actual(i) - expected(i))
      if (err > TOL) then
        write (*, '(A,1X,A,1X,I0,1X,ES12.4)') "mismatch", trim(label)//": index", i, err
        matches = .false.
        return
      end if
    end do
    write (*, '(A,1X,A)') "solution ok", trim(label)
  end function compare_solution

  subroutine copy_host_to_device(host, device_ptr)
    complex(C_DOUBLE_COMPLEX), intent(in), target :: host(:)
    type(C_PTR), intent(in) :: device_ptr
    integer(C_INT) :: status
    integer(C_SIZE_T) :: bytes

    bytes = int(storage_size(host(1))/8, kind=C_SIZE_T)*size(host, kind=C_SIZE_T)
    status = hipMemcpy(device_ptr, c_loc(host(1)), bytes, HIP_MEMCPY_HOST_TO_DEVICE)
    call check_hip_status(status, "hipMemcpy host_to_device")
  end subroutine copy_host_to_device

  subroutine copy_device_to_host(device_ptr, host)
    type(C_PTR), intent(in) :: device_ptr
    complex(C_DOUBLE_COMPLEX), intent(out), target :: host(:)
    integer(C_INT) :: status
    integer(C_SIZE_T) :: bytes

    bytes = int(storage_size(host(1))/8, kind=C_SIZE_T)*size(host, kind=C_SIZE_T)
    status = hipMemcpy(c_loc(host(1)), device_ptr, bytes, HIP_MEMCPY_DEVICE_TO_HOST)
    call check_hip_status(status, "hipMemcpy device_to_host")
  end subroutine copy_device_to_host

  subroutine check_status(status, where)
    integer(C_INT), intent(in) :: status
    character(*), intent(in) :: where
    if (status /= HIPSPARSE_STATUS_SUCCESS) then
      write (*, '(A,1X,A,1X,I0)') "hipSPARSE failure", trim(where)//":", status
      error stop
    end if
  end subroutine check_status

  subroutine check_hip_status(status, where)
    integer(C_INT), intent(in) :: status
    character(*), intent(in) :: where
    if (status /= 0_C_INT) then
      write (*, '(A,1X,A,1X,I0)') "HIP runtime failure", trim(where)//":", status
      error stop
    end if
  end subroutine check_hip_status

#endif

end program test_hipsparse_gpsv_allocations
