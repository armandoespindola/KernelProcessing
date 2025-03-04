! Ebru1: The choice of THRESHOLD value is somewhat subjective. It is not trivial to set it like the 20% of max value
! which may be OK for smaller scale studies but global scale needs a few trial&error to adjust this parameter for
! every iteration. Needs some more investigation..

! Ebru2: I find the preconditioner behave better after changing the order of smoothing and preconditioning in
! post-processing upon the suggestion by Ryan & Yanhua.
! However, I am still not convinced by Ryan's latest suggestion that preconditioner should be smoothed more than the
! gradients of other parameters that the preconditioner to be applied. I currently smooth the preconditioner and
! the other gradients in the same way.


module preconditioner_subs
  use mpi
  use global_var, only : max_all_all_cr, min_all_all_cr, CUSTOM_REAL, exit_mpi, &
       myrank
  use global_var,only : ATTENUATION_FLAG,KQMU_IDX,HESS_NAMES_GLOB,NHESS0
  implicit none

  contains

  subroutine get_sys_args(kernel_par,hess_file, output_dir, threshold_hess)
    character(len=*), intent(inout) :: output_dir, hess_file,kernel_par
    real(kind=CUSTOM_REAL), intent(inout) :: threshold_hess

    character(len=20) :: threshold_str

    call getarg(1, kernel_par)
    call getarg(2, hess_file)
    call getarg(3, output_dir)
    call getarg(4, threshold_str)

    if(hess_file == '' .or. output_dir == '' .or. threshold_str == '' &
         .or. kernel_par == '' ) then
      call exit_mpi("Usage: xinv_hess kernel_par hess_file outputdir  "// &
        "threshold")
    endif

    read(threshold_str, *) threshold_hess

    if (myrank == 0) then
       write(*, *) "Kernel Parfile: ", trim(kernel_par)
      write(*, *) "Hess kernel: ", trim(hess_file)
      write(*, *) "Preconditioner(output_dir): ", trim(output_dir)
      write(*, *) "Dampening threshold: ", threshold_hess
    endif

  end subroutine get_sys_args

  ! Take the inverse of Hessian with dampening (as preconditioner kernels)
  ! This method will apply different lambda for invidual model parameters
  ! but using the same threshold ratio regarding to its own max value
  subroutine prepare_invHess(hess, hess_names, threshold, invHess,invNames)

    real(CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: hess, invHess
    character(len=*), dimension(:),intent(in) :: hess_names
    character(len=*), dimension(:),intent(inout) :: invNames
    real(CUSTOM_REAL), intent(in) :: threshold

    real(kind=CUSTOM_REAL):: maxh_all, minh_all, dampening,ratio
    integer :: ik
    integer :: nhess

    nhess = int(size(hess_names))
    do ik = 1,nhess
       invNames(ik) = trim(hess_names(ik))
      if (myrank == 0) then
         write(*, '(/,A,A)') "====> Prepare invHess for parameter: ", &
              trim(invNames(ik))
      endif

      ! stats before
      call max_all_all_cr(maxval(hess(:, :, :, :, ik)), maxh_all)
      call min_all_all_cr(minval(hess(:, :, :, :, ik)), minh_all)

      if ( maxh_all < 1.e-18 ) then
        call exit_mpi("hess max value < 1.e-18")
      end if

      ! damping value
      dampening = maxh_all * threshold
      if (myrank == 0) then
        write(*, '(A, ES12.2, 5X, ES12.2, 5X, A, ES12.2)') &
          "Max and Min of hess: ", maxh_all, minh_all, &
          " || Condition Number:", maxh_all / minh_all
        write(*, '(A, E12.2, 5X, E12.2)') "Dampen threshold and value: ", threshold, dampening
      endif

      invHess(:, :, :, :, ik) = 1.0 / (hess(:, :, :, :, ik) + dampening)

      ! stats afterwards
      call max_all_all_cr(maxval(invHess(:, :, :, :, ik)), maxh_all)
      call min_all_all_cr(minval(invHess(:, :, :, :, ik)), minh_all)

      if (myrank == 0) then
        write(*, *) "Apply dampening and inverse..."
        write(*, '(A, ES12.2, 5X, ES12.2, 5X, A, ES12.2)') &
          "Max and Min of hess: ", maxh_all, minh_all, &
          " || Condition Number:", maxh_all / minh_all
      endif

   enddo

   if (myrank == 0) write(*,*) "|--------- Normalize Inv Hessian ---------|"
   call max_all_all_cr(maxval(invHess(:, :, :, :,:)), maxh_all)
   invHess = invHess / maxh_all

   if (myrank == 0) write(*,*) "|--------- Normalize Hessian ---------|"
   call max_all_all_cr(maxval(hess(:, :, :, :,:)), maxh_all)
   hess = hess / maxh_all
   
   call min_all_all_cr(minval(invHess(:, :, :, :,:)), minh_all)
   ratio = 1.0 + minh_all
   
   ! do ik = 1,nhess
   !    call max_all_all_cr(maxval(invHess(:, :, :, :,ik)), maxh_all)
   !    if (maxh_all < threshold**2) then
   !       invHess(:, :, :, :,ik) = invHess(:, :, :, :,ik) * (threshold**2) / (maxh_all)
   !    endif
   !    call max_all_all_cr(maxval(invHess(:, :, :, :,ik)), maxh_all)
   !    call min_all_all_cr(minval(invHess(:, :, :, :,ik)), minh_all)
   !    if (myrank == 0) then
   !       write(*, '(A, ES12.2, 5X, ES12.2, 5X, A, ES12.2)') &
   !            trim(invNames(ik)), maxh_all, minh_all, &
   !            " || Condition Number:", maxh_all / minh_all
   !    endif
   ! enddo

   if (ATTENUATION_FLAG) then
         !call max_all_all_cr(maxval(invHess(:, :, :, :,ik)), maxh_all)
         !call min_all_all_cr(minval(invHess(:, :, :, :,ik)), minh_all)
         write(*, *) "Multiparameter Inversion"
         invHess(:, :, :, :,KQMU_IDX) = invHess(:, :, :, :,KQMU_IDX) * 2.0e+3
         ! if ((ratio / (maxh_all + minh_all)) > 1000) then
         !    invHess(:, :, :, :,ik) = invHess(:, :, :, :,ik) * ratio * 1.0e-3 / (maxh_all + minh_all)
         ! endif
      endif
      call max_all_all_cr(maxval(invHess(:, :, :, :,:)), maxh_all)
      call min_all_all_cr(minval(invHess(:, :, :, :,:)), minh_all)
      invHess = invHess / maxh_all
      do ik = 1,nhess
         call max_all_all_cr(maxval(invHess(:, :, :, :,ik)), maxh_all)
         call min_all_all_cr(minval(invHess(:, :, :, :,ik)), minh_all)

      if (myrank == 0) then
         write(*, '(A, ES12.2, 5X, ES12.2, 5X, A, ES12.2)') &
              trim(invNames(ik)), maxh_all, minh_all, &
              " || Condition Number:", maxh_all / minh_all
      endif
   enddo
   
   call max_all_all_cr(maxval(invHess(:, :, :, :,:)), maxh_all)
   call min_all_all_cr(minval(invHess(:, :, :, :,:)), minh_all)

   if (myrank == 0) then
      write(*, *) "Apply dampening and inverse (after normalization)..."
      write(*, '(A, ES12.2, 5X, ES12.2, 5X, A, ES12.2)') &
           "Max and Min of hess: ", maxh_all, minh_all, &
           " || Condition Number:", maxh_all / minh_all
   endif
   
  end subroutine prepare_invHess

end module preconditioner_subs

program inverse_hessian
  use mpi
  use adios_read_mod
  use AdiosIO
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC, myrank, CUSTOM_REAL,KER_HESS_NAMES_GLOB
  use global_var, only : init_mpi,KERNEL_NAMES_GLOB,NKERNEL_GLOB,KERNEL_NAMES_GLOB,init_kernel_par
  use preconditioner_subs

  implicit none

  ! integer, parameter :: NHESS = 3
  ! character(len=500), parameter :: hess_names(NHESS) = &
  !   (/character(len=500) :: "hess_vp_kl_crust_mantle", &
  !                           "hess_vs_kl_crust_mantle", &
  !                           "hess_eta_kl_crust_mantle"/)
  ! real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NHESS) :: hess = 0.0, invHess = 0.0

  ! integer, parameter :: NPRECOND = 4
  ! character(len=500), parameter :: precond_names(NPRECOND) = &
  !   (/character(len=500) :: "precond_bulk_c_kl_crust_mantle", &
  !                           "precond_bulk_betav_kl_crust_mantle", &
  !                           "precond_bulk_betah_kl_crust_mantle", &
  !                           "precond_eta_kl_crust_mantle"/)
  ! real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NPRECOND) :: precond = 0.0

  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: hess,invHess
  character(len=500), dimension(:),allocatable :: precond_names
  character(len=500) :: hess_file, output_dir,kernel_parfile
  real(kind=CUSTOM_REAL) :: threshold_hess
  integer:: ier

  call init_mpi()
  
  if (myrank == 0) print*, "Prepare Preconditioners"

  call get_sys_args(kernel_parfile,hess_file, output_dir, threshold_hess)

  call init_kernel_par(kernel_parfile)

  allocate(hess(NGLLX, NGLLY, NGLLZ, NSPEC, NHESS0), &
       invHess(NGLLX, NGLLY, NGLLZ, NSPEC, NHESS0))

  allocate(precond_names(NHESS0))


  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
       "verbose=1", ier)

  
  call read_bp_file_real(hess_file, HESS_NAMES_GLOB, hess)

  call prepare_invHess(hess, HESS_NAMES_GLOB, threshold_hess, invHess,precond_names)

  call write_bp_file(invHess, precond_names, "KERNEL_GOURPS", trim(output_dir)//"invhess.bp")
  call write_bp_file(hess, precond_names, "KERNEL_GOURPS", trim(output_dir)//"hess.bp")
  
  if(myrank == 0) then
     write(*, '(/, A, A)') "output Inv Hess: ", trim(output_dir)//"invhess.bp"
     write(*, '(/, A, A)') "output Hess: ", trim(output_dir)//"hess.bp"
  endif

  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program inverse_hessian
