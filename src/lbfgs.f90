program main
  use mpi
  use adios_read_mod
  use lbfgs_subs
  use AdiosIO, only : write_bp_file, calculate_jacobian_matrix
  use global_var, only : init_mpi,NGLLX,NGLLY,NGLLZ,NPAR_GLOB,KERNEL_NAMES_GLOB,NSPEC,init_kernel_par,NKERNEL_GLOB
  implicit none

  !integer, parameter :: NKERNELS = 4
  !character(len=512), parameter :: kernel_names(NKERNELS) = &
  !      (/character(len=512) :: "bulk_c_kl_crust_mantle", &
  !                              "bulk_betav_kl_crust_mantle", &
  !                              "bulk_betah_kl_crust_mantle", &
  !                              "eta_kl_crust_mantle"/)

  ! Number of previous iteration used in L-BFGS
  integer :: niter
  character(len=512) :: input_path_file, solver_file, outputfn
  character(len=512), dimension(:), allocatable :: gradient_files, model_change_files

  ! jacobian related to the geometry of mesh
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: jacobian


  ! precond kernels (default = 1.0, no preconditioner applied)
  !!real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: precond = 1.0
  ! most recent gradient (this iteration's gradient)
  !!real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: gradient
  ! this iterations' search direction
  !!real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: direction

  ! precond kernels (default = 1.0, no preconditioner applied)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: precond
  ! most recent gradient (this iteration's gradient)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: gradient
  ! this iterations' search direction
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: direction

  ! According to Numerical Optimization, page 177:
  !     sk = x(k+1) - x(k), which is the model change
  !     yk = g(k+1) - g(k), which is the gradient change
  real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: yks, sks

  integer :: ier

  call init_mpi()
  call init_kernel_par()

  allocate(precond(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB),gradient(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB), &
       direction(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB))

  precond=1.0

  if(myrank == 0) print*, "|<---- Get System Args ---->|"
  call get_sys_args(input_path_file, solver_file, outputfn)

  if(myrank == 0) print*, "|<---- Parse Input Path File ---->|"
  call parse_input_path(input_path_file, niter, gradient_files, model_change_files)

  allocate(sks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB, niter))
  allocate(yks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB, niter))

  call read_all_bp_files(niter, NKERNEL_GLOB, gradient_files, model_change_files, &
                         KERNEL_NAMES_GLOB, gradient, yks, sks)

  if(myrank == 0) print*, "|<---- Calculate Jacobian ---->|"
  call calculate_jacobian_matrix(solver_file, jacobian)

  if(myrank == 0) print*, "|<---- L-BFGS Direction ---->|"
  call calculate_LBFGS_direction(niter, NKERNEL_GLOB, jacobian, gradient, precond, yks, sks, direction)

  call write_bp_file(direction, KERNEL_NAMES_GLOB, "KERNELS_GROUP", outputfn)
  if(myrank == 0) print*, "LBFGS direction saved: ", trim(outputfn)

  call adios_finalize(myrank, ier)

  deallocate(precond,gradient,direction)
  
  call MPI_finalize(ier)

end program main
