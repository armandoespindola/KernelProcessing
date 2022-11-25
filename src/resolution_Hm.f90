program main
  use mpi
  use adios_read_mod
  use lbfgs_subs
  use AdiosIO, only : write_bp_file, calculate_jacobian_matrix
  use global_var, only : init_mpi,NGLLX,NGLLY,NGLLZ,NPAR_GLOB,KERNEL_NAMES_GLOB,NSPEC,init_kernel_par,NKERNEL_GLOB,HESS_NAMES_GLOB,max_all_all_cr,min_all_all_cr
  implicit none

  !integer, parameter :: NKERNELS = 4
  !character(len=512), parameter :: kernel_names(NKERNELS) = &
  !      (/character(len=512) :: "bulk_c_kl_crust_mantle", &
  !                              "bulk_betav_kl_crust_mantle", &
  !                              "bulk_betah_kl_crust_mantle", &
  !                              "eta_kl_crust_mantle"/)

  ! Number of previous iteration used in L-BFGS
  integer :: niter
  character(len=512) :: input_path_file, solver_file, outputdir,kernel_parfile,precond_file,perturb_file,output_file
  character(len=512) :: grad_m0_file, grad_dm_file
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
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: precond,diagH
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: gradient_dm, gradient_m0,test_m
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: Bm

  ! According to Numerical Optimization, page 177:
  !     sk = x(k+1) - x(k), which is the model change
  !     yk = g(k+1) - g(k), which is the gradient change
  real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: yks, sks
  real(kind=CUSTOM_REAL) :: max_direction,min_direction

  integer :: ier,iker

  call init_mpi()

  if(myrank == 0) print*, "|<---- Get System Args ---->|"
  call get_sys_args_bm(kernel_parfile,input_path_file, solver_file, grad_dm_file, &
       precond_file, outputdir)

  call init_kernel_par(kernel_parfile)

  allocate(precond(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB), &
       gradient_dm(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB), &
       Bm(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB), &
       test_m(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB))


  if(myrank == 0) print*, "|<---- Parse Input Path File ---->|"
  call parse_input_path(input_path_file, niter, gradient_files, model_change_files)

  allocate(sks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB, niter))
  allocate(yks(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB, niter))

  call read_all_bp_files(niter, NKERNEL_GLOB, gradient_files, model_change_files, &
                         KERNEL_NAMES_GLOB, test_m, yks, sks)

  test_m = 0.0d0
  
  if(myrank == 0) print*, "|<---- Calculate Jacobian ---->|"
  call calculate_jacobian_matrix(solver_file, jacobian)
  jacobian = 1.0d0

  if (trim(precond_file) == '') then
     precond = 1.0
     if (myrank == 0) write(*, *) "|<----- Preconditioner H=I ----->| "
  else
     if (myrank == 0) write(*, *) "|<----- Reading Preconditioner ----->| ", &
          trim(precond_file)
     call read_bp_file_real(precond_file, HESS_NAMES_GLOB, precond)
  endif

  
  if(myrank == 0) print*, "|<---- Reading dm ---->|"
  call read_bp_file_real(grad_dm_file, KERNEL_NAMES_GLOB, gradient_dm)

  test_m = gradient_dm !- gradient_m0

  
  if(myrank == 0) print*, "|<---- L-BFGS Bm ---->|"
  call calculate_LBFGS_Bm(niter, NKERNEL_GLOB, jacobian, test_m, precond, yks, sks, Bm)

  ! if(myrank == 0) print*, "|<---- L-BFGS H^-1(Bdm) ---->|"
  ! call calculate_LBFGS_direction(niter, NKERNEL_GLOB, jacobian, test_m, precond, yks, sks, Bm)

  if(myrank == 0) print*, "|<---- L-BFGS Bm (Stats) ---->|"
  do iker = 1,NKERNEL_GLOB
     call max_all_all_cr(maxval(Bm(:, :, :, :,iker)), max_direction)
     call min_all_all_cr(minval(Bm(:, :, :, :,iker)), min_direction)
     if (myrank == 0) then
        write(*,*) " Bm Stats "
        write(*, '(A, ES12.2, 5X, ES12.2)') &                 
             trim(KERNEL_NAMES_GLOB(iker))//"(max,min)", & 
             max_direction, min_direction                
     endif
  enddo


  ! output_file = trim(outputdir)//'/Bdm.bp'
  ! call write_bp_file(test_m, KERNEL_NAMES_GLOB, "KERNELS_GROUP", output_file )
  ! if(myrank == 0) print*, "Resolution Bdm saved: ", trim(output_file)

  output_file = trim(outputdir)//'/Bdm.bp'
  call write_bp_file(Bm, KERNEL_NAMES_GLOB, "KERNELS_GROUP", output_file)
  if(myrank == 0) print*, "Resolution LBFGS (Bdm) saved: ", trim(output_file)

  call adios_finalize(myrank, ier)

  deallocate(precond,Bm,test_m)
  
  call MPI_finalize(ier)

end program main
