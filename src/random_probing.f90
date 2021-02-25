! -----------------------------------------------------------------------
! Random probing for L-BFGS
!
! It will feed a vector, filled with random values from normal distribution,
! feed them into the the L-BFGS loop
!
! Input:
!   1) L-BFGS loop path file, including the past gradients and module changes
!   2) solver file, to caculate the jacobian
!   3) precondition bp file: the precondition kernels are inverse
!         Hessian (approximate). You need to inverse it.
! -----------------------------------------------------------------------


module random_probe_mod

  use mpi

  use global_var, only : CUSTOM_REAL, myrank, exit_mpi, NGLLX, NGLLY, NGLLZ, NSPEC
  use AdiosIO, only : read_bp_file_real

  implicit none

  contains

  subroutine get_sys_args(lbfgs_path_file, solver_file, precond_file, output_file_prefix)

    character(len=*), intent(out) :: lbfgs_path_file, solver_file, precond_file
    character(len=*), intent(out) :: output_file_prefix

    call getarg(1, lbfgs_path_file)
    call getarg(2, solver_file)
    call getarg(3, precond_file)
    call getarg(4, output_file_prefix)

    if(trim(lbfgs_path_file) == "" .or. trim(solver_file) == "" .or. &
      trim(precond_file) == "" .or. trim(output_file_prefix) == "") then
      call exit_mpi("Usage: ./xrand_probe lbfgs_path_file solver_file precond_file output_file_prefix")
    endif

    if(myrank == 0) then
      print*, "Input L-BFGS path: ", trim(lbfgs_path_file)
      print*, "Solver file: ", trim(solver_file)
      print*, "Precond file: ", trim(precond_file)
      print*, "Output file prefix: ", trim(output_file_prefix)
    endif
  end subroutine get_sys_args


  subroutine read_lbfgs_history(niter, nkernels, kernel_names, lbfgs_path_file, yks, sks)
    use lbfgs_subs, only : parse_input_path, read_all_bp_files

    integer, intent(out) :: niter
    integer, intent(in) :: nkernels
    character(len=*), dimension(:) :: kernel_names

    character(len=*), intent(in) :: lbfgs_path_file

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: gradient
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: yks, sks

    character(len=512), dimension(:), allocatable :: gradient_files, model_change_files

    allocate(sks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS, niter))
    allocate(yks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS, niter))

    call parse_input_path(lbfgs_path_file, niter, gradient_files, model_change_files)

    call read_all_bp_files(niter, nkernels, gradient_files, model_change_files, &
                           kernel_names, gradient, yks, sks)

  end subroutine read_lbfgs_history


  subroutine prepare_preconditioner(nkernels, nhess, hess_names, precond_file, precond)

    use precond_vp_vs_kernels_subs, only : prepare_hessian

    integer, intent(in) :: nkernels, nhess
    character(len=*), dimension(:) :: hess_names
    character(len=*), intent(in) :: precond_file
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nkernels), intent(out) :: precond

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, nhess) :: hess, invHess
    real(kind=CUSTOM_REAL) :: threshold = 1.0e-3

    hess = 0.0
    invHess = 0.0

    call read_bp_file_real(precond_file, hess_names, hess)
    call prepare_hessian(hess, hess_names, threshold, invHess)

    ! bulk_c -- Vp precondition
    precond(:, :, :, :, 1) = invHess(:, :, :, :, 1)

    ! bulk_betav and bulkd_betah -- Vs preconditioner
    precond(:, :, :, :, 2) = invHess(:, :, :, :, 2)
    precond(:, :, :, :, 3) = invHess(:, :, :, :, 2)

    ! eta -- eta preconditioner
    precond(:, :, :, :, 4) = invHess(:, :, :, :, 3)

  end subroutine prepare_preconditioner


  subroutine random_fill(nkernel, r)

    use random, only : random_normal

    integer, intent(in) :: nkernel
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(out) :: r

    integer :: iker, ispec, i, j, k

    do iker = 1, nkernel
      do ispec = 1, NSPEC
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i =1, NGLLX
              r(i, j, k, ispec, iker) = random_normal()
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine random_fill


end module random_probe_mod

program main
  use mpi
  use adios_read_mod
  use AdiosIO, only : write_bp_file, calculate_jacobian_matrix
  use global_var, only : CUSTOM_REAL, myrank, init_mpi, exit_mpi
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC

  use random_probe_mod, only : get_sys_args, read_lbfgs_history, &
    prepare_preconditioner, random_fill
  use lbfgs_subs, only : calculate_LBFGS_direction

  implicit none

  integer, parameter :: NKERNELS = 4
  character(len=512), parameter :: kernel_names(NKERNELS) = &
        (/character(len=512) :: "bulk_c_kl_crust_mantle", &
                                "bulk_betav_kl_crust_mantle", &
                                "bulk_betah_kl_crust_mantle", &
                                "eta_kl_crust_mantle"/)

 integer, parameter :: NHESS = 3
  character(len=512), parameter :: hess_names(NHESS) = &
    (/character(len=512) :: "hess_vp_kl_crust_mantle", "hess_vs_kl_crust_mantle", &
                            "hess_eta_kl_crust_mantle"/)

  integer :: niter
  character(len=512) :: lbfgs_path_file, solver_file, precond_file
  character(len=512) :: output_file_prefix, output_file

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: jacobian

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS) :: precond
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS) :: r1, r2

  real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: yks, sks

  integer :: ier

  call init_mpi()

  if(myrank == 0) print*, "|<---- Get System Args ---->|"
  call get_sys_args(lbfgs_path_file, solver_file, precond_file, output_file_prefix)

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  if(myrank == 0) print*, "|<---- Get System Args ---->|"
  call read_lbfgs_history(niter, nkernels, kernel_names, lbfgs_path_file, yks, sks)

  if(myrank == 0) print*, "|<---- Calculate Jacobian ---->|"
  call calculate_jacobian_matrix(solver_file, jacobian)

  if(myrank == 0) print*, "|<---- Calculate Jacobian ---->|"
  call prepare_preconditioner(NKERNELS, NHESS, hess_names, precond_file, precond)

  call random_fill(NKERNELS, r1)
  call random_fill(NKERNELS, r2)

  ! --- l-bfgs random probing ---
  if(myrank == 0) print*, "|<---- Random Probing in L-BFGS loop ---->|"
  call calculate_LBFGS_direction(niter, NKERNELS, jacobian, r1, yks, sks, r2)
  r2 = r2 - r1

  ! --- write ---
  output_file = trim(output_file_prefix) // "_R.bp"
  call write_bp_file(r1, kernel_names, "KERNELS_GROUP", output_file)

  output_file = trim(output_file_prefix) // "_r.bp"
  call write_bp_file(r2, kernel_names, "KERNELS_GROUP", output_file)

  deallocate(yks)
  deallocate(sks)

  call adios_finalize(myrank, ier)
  call mpi_finalize(ier)

end program main
