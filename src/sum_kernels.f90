module sum_kernels_subs

  use global_var, only : CUSTOM_REAL, myrank, exit_mpi
  implicit none

  contains

  subroutine read_event_file(eventfile, nevent, kernel_list, weight_list)
    character(len=*), intent(in) :: eventfile
    integer, intent(inout) :: nevent
    character(len=*), dimension(:), allocatable, intent(inout) :: kernel_list
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: weight_list

    ! local variables
    character(len=500) :: kernel_file
    real(kind=CUSTOM_REAL) :: weight
    integer :: ios, i
    integer, parameter :: fh = 1001

    open(unit=fh, file=trim(eventfile), status='old', iostat=ios)
    if ( ios /= 0 ) then
      print*, 'ERROR OPENING', trim(eventfile)
      stop
    end if

    if(myrank == 0) write(*, *) "Reading events and weights from: ", trim(eventfile)
    read(fh, *) nevent
    if(myrank == 0) write(*, *) "Number of events: ", nevent

    allocate(kernel_list(nevent))
    allocate(weight_list(nevent))

    do i=1, nevent
      read(fh, *) weight
      read(fh, '(A)') kernel_file
      if (myrank == 0) write(*, '(A, I5, A, I5, A, A, A, ES16.8)') &
        "[", i, "/", nevent, "]", trim(kernel_file), ": ", weight
      kernel_list(i) = trim(kernel_file)
      weight_list(i) = weight
    enddo
    close(fh)
  end subroutine read_event_file

  subroutine get_sys_args(eventfile, outputfn)
    character(len=*), intent(inout) :: eventfile, outputfn

    call getarg(1, eventfile)
    call getarg(2, outputfn)

    if(trim(eventfile) == '' .or. trim(outputfn) == '') then
      call exit_mpi("Usage: xsum_kernels eventfile outputfn")
    endif

    if(myrank == 0) then
      write(*, *) "Event file(in): ", trim(eventfile)
      write(*, *) "Output file(out, kernel sums): ", trim(outputfn)
    endif
  end subroutine

end module sum_kernels_subs


program sum_kernels
  use mpi
  use adios_read_mod
  use global_var, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank,init_kernel_par,NPAR_GLOB,NKERNEL_GLOB
  use global_var, only : init_mpi, exit_mpi,KERNEL_NAMES_GLOB,MODEL_NAMES_GLOB
  use AdiosIO
  use sum_kernels_subs

  implicit none

  integer:: nevent, ievent, ier
  character(len=500) :: eventfile, outputfn, kernel_file

  character(len=500), dimension(:), allocatable :: kernel_list
  real(kind=CUSTOM_REAL):: weight
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: weight_list

  !integer, parameter :: NKERNELS = 5    !bulk_betah, bulk_betav, bulk_c, eta, rho, hess
  !integer, parameter :: hess_idx = 5
  !character(len=500), parameter :: kernel_names(NKERNELS) = &
  !  (/character(len=150) :: "bulk_betah_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
  !                          "bulk_c_kl_crust_mantle",     "eta_kl_crust_mantle",        &
  !                          "hess_kl_crust_mantle"/)

  ! integer, parameter :: NKERNELS = 7    !bulk_betah, bulk_betav, bulk_c, eta, rho, hess
  ! integer, parameter :: hess_idx = 5
  ! character(len=500), parameter :: kernel_names(NKERNELS) = &
  !   (/character(len=500) :: "bulk_betah_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
  !                           "bulk_c_kl_crust_mantle", "eta_kl_crust_mantle",        &
  !                           "hess_kl_crust_mantle", "hess_kappa_kl_crust_mantle", &
  !                           "hess_mu_kl_crust_mantle"/)

  ! real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: total_kernel
  ! real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: kernels


  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: total_kernel
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: kernels

  call init_mpi()

  call init_kernel_par()

  allocate(total_kernel(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB),kernels(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB))

!  if (kernel_names(hess_idx) /= "hess_kl_crust_mantle") call exit_mpi("Incorrect hess_idx")

  call get_sys_args(eventfile, outputfn)
  call read_event_file(eventfile, nevent, kernel_list, weight_list)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  total_kernel = 0.
  do ievent=1, nevent
    kernel_file = kernel_list(ievent)
    weight = weight_list(ievent)

    if (myrank==0) write(*,*) 'Reading in kernel for [', ievent, &
        "/", nevent, "]: ", trim(kernel_file), " | weight: ", weight

    call read_bp_file_real(kernel_file, KERNEL_NAMES_GLOB, kernels)

    ! only keep the abs value of the hess kernel
    !kernels(:, :, :, :, hess_idx) = abs(kernels(:, :, :, :, hess_idx))
    !kernels(:, :, :, :, hess_idx:(hess_idx+2)) = abs(kernels(:, :, :, :, hess_idx:(hess_idx+2)))

    total_kernel = total_kernel + kernels * weight
  end do

  call write_bp_file(total_kernel, KERNEL_NAMES_GLOB, "KERNEL_GROUPS", outputfn)

  if (myrank==0) print*, 'Done summing all the kernels'
  if (myrank==0) print*, "output file: ", trim(outputfn)

  call adios_finalize(myrank, ier)
  deallocate(total_kernel,kernels)
  call MPI_FINALIZE(ier)

end program sum_kernels
