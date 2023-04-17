! add_model_globe_tiso
!
! this program can be used to update TRANSVERSE ISOTROPIC model files
! based on smoothed event kernels.
! the kernels are given for tranverse isotropic parameters (bulk_c,bulk_betav,bulk_betah,eta).
!
! the algorithm uses a steepest descent method with a step length
! determined by the given maximum update percentage.
!
! input:
!    - step_fac : step length to update the models, f.e. 0.03 for plusminus 3%
!
! setup:
!
!- INPUT_MODEL/  contains:
!       proc000***_reg1_vsv.bin &
!       proc000***_reg1_vsh.bin &
!       proc000***_reg1_vpv.bin &
!       proc000***_reg1_vph.bin &
!       proc000***_reg1_eta.bin &
!       proc000***_reg1_rho.bin
!
!- INPUT_GRADIENT/ contains:
!       proc000***_reg1_bulk_c_kernel_smooth.bin &
!       proc000***_reg1_bulk_betav_kernel_smooth.bin &
!       proc000***_reg1_bulk_betah_kernel_smooth.bin &
!       proc000***_reg1_eta_kernel_smooth.bin
!
!- topo/ contains:
!       proc000***_reg1_solver_data_1.bin
!
! new models are stored in
!- OUTPUT_MODEL/ as
!   proc000***_reg1_vpv_new.bin &
!   proc000***_reg1_vph_new.bin &
!   proc000***_reg1_vsv_new.bin &
!   proc000***_reg1_vsh_new.bin &
!   proc000***_reg1_eta_new.bin &
!   proc000***_reg1_rho_new.bin
!
! USAGE: ./add_model_globe_tiso 0.3

! 1st Version by Hejun Zhu
! 2nd Version with Adios implementation: Ebru & Matthieu, August 2013
! 3rd Version by Wenjie Lei, June 2018
! Princeton

module model_update_par

  use mpi
  use global_var, only : myrank, nprocs, NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB, CUSTOM_REAL
  use global_var, only : max_all_all_cr, min_all_all_cr, init_mpi, exit_mpi
  use global_var, only : KERNEL_NAMES_GLOB,MODEL_NAMES_GLOB,MODEL_PERTURB_NAMES_GLOB,NPAR_GLOB,NKERNEL_GLOB,NMODEL0,MODEL0_NAMES
  use global_var, only : init_kernel_par,NMODEL_TOTAL,MODEL_NAMES_TOTAL
  use global_var, only : VPV_IDX,VPH_IDX,VSV_IDX,VSH_IDX,RHO_IDX,ETA_IDX
  use global_var, only : KVPV_IDX,KVPH_IDX,KVSV_IDX,KVSH_IDX,KRHO_IDX,KETA_IDX,ATTENUATION_FLAG
  use AdiosIO

  implicit none

  logical, parameter :: use_depth_maximum = .false.  ! false
  ! ======================================================
  ! density scaling factor with shear perturbations
  ! see e.g. Montagner & Anderson (1989), Panning & Romanowicz (2006)
  real(kind=CUSTOM_REAL), parameter :: RHO_SCALING   = 0.33_CUSTOM_REAL
  ! constraint on eta model
  real(kind=CUSTOM_REAL), parameter :: LIMIT_ETA_MIN = 0.5_CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: LIMIT_ETA_MAX = 1.5_CUSTOM_REAL

  ! ! ======================================================
  ! ! MODELS
  ! integer, parameter :: NMODELS = 6
  ! character(len=500), dimension(NMODELS), parameter :: model_names = &
  !   (/character(len=500) :: "reg1/vpv", "reg1/vph", "reg1/vsv", &
  !                           "reg1/vsh", "reg1/eta", "reg1/rho"/)
  ! integer, parameter :: vpv_idx=1, vph_idx=2, vsv_idx=3, vsh_idx=4, &
  !                       eta_idx=5, rho_idx=6
  ! character(len=500), dimension(NMODELS), parameter :: MODEL_PERTURB_NAMES_GLOB = &
  !   (/character(len=150) :: "reg1/dvpvvpv","reg1/dvphvph","reg1/dvsvvsv", &
  !                           "reg1/dvshvsh","reg1/detaeta","reg1/drhorho"/)
  ! ! transverse isotropic model files
  ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NMODELS) :: models = 0.0
  ! real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NMODELS) :: models_new = 0.0

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: models,models_fixed
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: models_new

  ! ! ======================================================
  ! ! KERNELS
  ! integer, parameter :: NKERNELS = 4
  ! character(len=500), dimension(NKERNELS), parameter :: kernel_names = &
  !   (/character(len=150) :: "bulk_c_kl_crust_mantle", &
  !                           "bulk_betav_kl_crust_mantle", &
  !                           "bulk_betah_kl_crust_mantle", &
  !                           "eta_kl_crust_mantle"/)
  ! integer, parameter :: bulk_c_kl_idx = 1, betav_kl_idx = 2, betah_kl_idx = 3, &
  !      eta_kl_idx = 4

  !   real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: kernels = 0.0
  ! ! model updates
  !   real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNELS) :: dmodels = 0.0
    
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:) , allocatable:: kernels
  ! model updates
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: dmodels

  ! ======================================================
  ! MESH information
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: x_glob, y_glob, z_glob
  integer, dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: ibool
  integer, dimension(NSPEC) :: idoubling
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: jacobian
  integer :: npar_m,nkernel

contains


  subroutine init_module_par(kernel_parfile)
    character(len=*), intent(inout) :: kernel_parfile
    integer :: ier,n_newmodel

    ! Initialize parameter Kernel Processing
    call init_kernel_par(kernel_parfile)
    
    ! Allocates memory for variables in model_par
    allocate(models(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB), &
         models_new(NGLLX,NGLLY,NGLLZ,NSPEC,NMODEL_TOTAL), &
         kernels(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB), &
         dmodels(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB),stat=ier)

    models=0.0d0
    models_new=0.0d0
    kernels=0.0d0
    dmodels=0.0d0

    if (ier /= 0) stop " Error Allocating parameter for Kernel Processing: Update Model "
    
    if (NMODEL0 .gt. 0) then
       allocate(models_fixed(NGLLX,NGLLY,NGLLZ,NSPEC,NMODEL0),stat=ier)
       if (ier /= 0) stop " Error Allocating parameter for Kernel Processing: Model Fixed "
    endif
    
    
    

  end subroutine init_module_par


  subroutine clean_module_par()

    ! Clean memory arrays
    deallocate(models,models_new,kernels,dmodels)
    
    end subroutine clean_module_par
  

  subroutine get_sys_args(kernel_parfile,step_fac, input_model_file, input_solver_file, &
                          input_kernel_file, outputdir)
    ! input args
    real(kind=CUSTOM_REAL), intent(inout) :: step_fac
    character(len=*), intent(inout) :: input_model_file, input_solver_file, &
                                       input_kernel_file, outputdir,kernel_parfile

    ! local
    character(len=150) :: s_step_fac

    call getarg(1, kernel_parfile)
    call getarg(2, s_step_fac)
    call getarg(3, input_model_file)
    call getarg(4, input_solver_file)
    call getarg(5, input_kernel_file)
    call getarg(6, outputdir)

    if (trim(s_step_fac) == '' .or. trim(input_model_file) == '' &
        .or. trim(input_kernel_file) == ''.or. &
         trim(input_solver_file) == '' .or. trim(outputdir) == '' .or. trim(kernel_parfile) == '') then
        call exit_MPI('Usage: add model_globe_tiso step_factor input_model input_solver input_kernel outputdir')
    endif

    read(s_step_fac, *) step_fac

    if( abs(step_fac) < 1.e-10) then
      if(myrank == 0) print *, 'error: step factor ', step_fac
      call exit_MPI('error step factor')
    endif

    if (myrank == 0) then
      print*,'Model update for vsv,vsh,vpv,vph,eta,rho'
      print*, "System Args: "
      write(*, '(A, ES16.6)'),' step_fac = ', step_fac
      print*,' input model file  :', trim(input_model_file)
      print*,' input solver file :', trim(input_solver_file)
      print*,' input kernel file :', trim(input_kernel_file)
      print*,' outputdir  :', trim(outputdir)
      print*
    endif
  end subroutine get_sys_args

  subroutine read_solver_data(input_solver_file)
    use AdiosIO, only : calculate_jacobian_matrix
    use global_var, only : Parallel_ComputeIntegral, sum_all_all_cr

    character(len=*), intent(in) :: input_solver_file

    integer :: i
    real(kind=CUSTOM_REAL) :: integl, norml
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: tmp_array

    call read_bp_file_int(input_solver_file, "reg1/idoubling", idoubling)
    call read_bp_file_int(input_solver_file, "reg1/ibool", ibool)
    call read_bp_file_real(input_solver_file, "reg1/x_global", x_glob)
    call read_bp_file_real(input_solver_file, "reg1/y_global", y_glob)
    call read_bp_file_real(input_solver_file, "reg1/z_global", z_glob)

    
    ! calculate jacobian matrix
    call calculate_jacobian_matrix(input_solver_file, jacobian)

    ! integral
    if(myrank == 0) print*, "Integral of Kernels: "
    do i = 1, nkernel
      call Parallel_ComputeIntegral(kernels(:, :, :, :, i), jacobian, integl)
      if(myrank == 0) write(*, '(4X,A30,A,ES16.8)') trim(KERNEL_NAMES_GLOB(i)), ": ", integl
    enddo
    tmp_array = 1.0
    call Parallel_ComputeIntegral(tmp_array, jacobian, integl)
    if(myrank == 0) write(*, '(4X,A32,ES16.8)') "Total Volume: ", integl

    ! norm
    if(myrank == 0) print*, "Norm of Kernels:"
    do i = 1, nkernel
      call sum_all_all_cr(sum(kernels(:, :, :, :, i) * kernels(:, :, :, :, i)), norml)
      norml = sqrt(norml)
      if(myrank == 0) write(*, '(A30, A, ES16.8)') trim(KERNEL_NAMES_GLOB(i)), ": ", norml
    enddo

  end subroutine read_solver_data

  subroutine get_model_change(step_fac)
    use global_var, only : R_EARTH_KM,NPAR_GLOB
    implicit none

    real(kind=CUSTOM_REAL), intent(inout) :: step_fac
    ! local parameters
    ! ------------------------------------------------------------------------
    ! sets maximum update in this depth range
    ! normalized radii
    double precision, parameter :: R_top = (R_EARTH_KM - 50.0 ) / R_EARTH_KM
    double precision, parameter :: R_bottom = (R_EARTH_KM - 600.0 ) / R_EARTH_KM
    real(kind=CUSTOM_REAL) :: r, vmax_depth, betav,global_vmax_all

    integer :: i, j, k, ispec, iglob

    real(kind=CUSTOM_REAL), dimension(:),allocatable:: vmax,global_vmax

    allocate(vmax(nkernel),global_vmax(nkernel))
    
    ! vmax = 0.0
    ! global_vmax = 0.0
    ! vmax_depth = 0.0
    ! ! gradient in negative direction for steepest descent
    ! ! determines maximum kernel betav value within given radius
    ! if( use_depth_maximum ) then
    !   if(myrank == 0) write(*, *) 'Using depth maximum between 50km and 600km'
    !   do ispec = 1, NSPEC
    !     do k = 1, NGLLZ
    !       do j = 1, NGLLY
    !         do i = 1, NGLLX
    !           ! get radius of point
    !           iglob = ibool(i,j,k,ispec)
    !           r = sqrt(x_glob(iglob)*x_glob(iglob) + y_glob(iglob)*y_glob(iglob) + &
    !                    z_glob(iglob)*z_glob(iglob))
    !           ! stores maximum kernel betav value in this depth slice
    !           ! since betav is most likely dominating
    !           if( r < R_top .and. r > R_bottom ) then
    !             ! kernel betav value
    !             betav = abs( kernels(i, j, k, ispec, betav_kl_idx) )
    !             if( vmax < betav ) then
    !               vmax = betav
    !               vmax_depth = r
    !             endif
    !           endif
    !         enddo
    !       enddo
    !     enddo
    !   enddo
    ! else
    !   ! maximum of all processes stored in max_vsv
    !   vmax = maxval(abs(kernels(:, :, :, :, betav_kl_idx)))
    !   if(myrank == 0) then
    !     write(*, *) 'Using vsv(all depth) as maximum'
    !   endif
    ! endif


    do i=1,nkernel
       vmax(i) = maxval(abs(kernels(:, :, :, :, i)))
    end do

       
    if(myrank == 0) print*, "Initial Model Update (before norm): "
    call stats_value_range(kernels, KERNEL_NAMES_GLOB)

    do i=1,nkernel
       call max_all_all_cr(vmax(i), global_vmax(i))
    enddo

    global_vmax_all = maxval(global_vmax)
    
    if(myrank == 0) then
      !print*, "--- Normalization ---"
      !write(*, '(A, ES16.8)') 'Using maximum to normalize: ', global_vmax_all
      write(*, '(A, F16.8)') 'Step length:                ', step_fac
    endif

!    dmodels = step_fac * (kernels / global_vmax_all)
    dmodels = step_fac * kernels
    
    if(myrank == 0) print*, "Scaled Model Update: "
    call stats_value_range(dmodels, KERNEL_NAMES_GLOB)


    deallocate(vmax,global_vmax)

  end subroutine get_model_change

  
  subroutine model_fixed()
    integer :: i
    do i=NPAR_GLOB+1,NMODEL_TOTAL
       models_new(:,:,:,:,i) = models_fixed(:,:,:,:,i - NPAR_GLOB)
    enddo
  end subroutine model_fixed

  subroutine update_model_tiso()
    
    if ((VPV_IDX .gt. 0) .and. (VPH_IDX .gt. 0)) call update_vp_tiso()
    if ((VSH_IDX .gt. 0) .and. (VSV_IDX .gt. 0)) call update_vs_tiso()
    if (ETA_IDX .gt. 0) call update_eta()
    if (RHO_IDX .gt. 0) call update_density()

    if (ATTENUATION_FLAG) call update_qmu()

    if(myrank == 0) print *, "Old Model:"
    call stats_value_range(models, MODEL_NAMES_GLOB)

    if (NMODEL0 .gt. 0) then
       call model_fixed()
       if(myrank == 0) print *, "New Model:"
       call stats_value_range(models_new, MODEL_NAMES_TOTAL)
    else
       if(myrank == 0) print *, "New Model:"
       call stats_value_range(models_new, MODEL_NAMES_TOTAL)
    endif
    
    
  end subroutine update_model_tiso


  subroutine update_density()
    use global_var , only : FOUR_THIRDS, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220
    ! model update:
    ! transverse isotropic update only in layer Moho to 220
    ! (where SPECFEM3D_GLOBE considers TISO)
    ! everywhere else uses an isotropic update
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
      model_dbetah, model_dbetav

    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL) :: betav0,betah0,rho0
    real(kind=CUSTOM_REAL) :: betav1,betah1,rho1
    real(kind=CUSTOM_REAL) :: betaiso0,betaiso1
    real(kind=CUSTOM_REAL) :: dbetaiso


    model_dbetav = dmodels(:, :, :, :, KVSV_IDX)
    model_dbetah = dmodels(:, :, :, :, KVSH_IDX)

    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              
            ! initial model values
            betav0  = models(i,j,k,ispec,VSV_IDX)
            betah0  = models(i,j,k,ispec,VSH_IDX)
            rho0  = models(i,j,k,ispec,RHO_IDX)
            rho1    = 0._CUSTOM_REAL

            ! do not use transverse isotropy except if element is between d220 and Moho
            if(.not. ( idoubling(ispec)==IFLAG_670_220 .or. idoubling(ispec)==IFLAG_220_80 &
                       .or. idoubling(ispec)==IFLAG_80_MOHO) ) then

              ! shear values
              ! isotropic kernel K_beta = K_betav + K_betah
              ! with same scaling step_length the model update dbeta_iso = dbetav + dbetah
              ! note:
              !   this step length can be twice as big as that given by the input
              dbetaiso = model_dbetav(i,j,k,ispec) + model_dbetah(i,j,k,ispec)
              ! note: betah is probably not really used in isotropic layers
              !         (see SPECFEM3D_GLOBE/get_model.f90)

              ! density: uses scaling relation with isotropic shear perturbations
              !               dln rho = RHO_SCALING * dln betaiso
              rho1 = rho0 * exp( RHO_SCALING * dbetaiso )
              
            else
              
              ! shear values
              betav1 = betav0 * exp( model_dbetav(i,j,k,ispec) )
              betah1 = betah0 * exp( model_dbetah(i,j,k,ispec) )

              ! density: uses scaling relation with Voigt average of shear perturbations
              betaiso0 = sqrt(  ( 2.0 * betav0**2 + betah0**2 ) / 3.0 )
              betaiso1 = sqrt(  ( 2.0 * betav1**2 + betah1**2 ) / 3.0 )
              dbetaiso = log( betaiso1 / betaiso0 )
              rho1 = rho0 * exp( RHO_SCALING * dbetaiso )
              
            endif
            ! stores new model values
            models_new(i,j,k,ispec,RHO_IDX) = rho1
          enddo
        enddo
      enddo
    enddo
  end subroutine update_density

  subroutine update_eta()
    use global_var , only : FOUR_THIRDS, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220
    ! model update:
    ! transverse isotropic update only in layer Moho to 220
    ! (where SPECFEM3D_GLOBE considers TISO)
    ! everywhere else uses an isotropic update
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
      model_deta

    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL) :: eta0
    real(kind=CUSTOM_REAL) :: eta1


    model_deta   = dmodels(:, :, :, :, KETA_IDX)

    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! initial model values
            eta0    = models(i,j,k,ispec,ETA_IDX)
            eta1    = 0._CUSTOM_REAL

            ! do not use transverse isotropy except if element is between d220 and Moho
            if(.not. ( idoubling(ispec)==IFLAG_670_220 .or. idoubling(ispec)==IFLAG_220_80 &
                       .or. idoubling(ispec)==IFLAG_80_MOHO) ) then
              ! isotropic model update
              ! no eta perturbation, since eta = 1 in isotropic media
              eta1 = eta0
            else
              ! transverse isotropic model update
              ! eta value : limits updated values for eta range constraint
              eta1 = eta0 * exp( model_deta(i,j,k,ispec) )

              if( eta1 < LIMIT_ETA_MIN ) eta1 = LIMIT_ETA_MIN
              if( eta1 > LIMIT_ETA_MAX ) eta1 = LIMIT_ETA_MAX
            endif
            ! stores new model values
            models_new(i,j,k,ispec,ETA_IDX) = eta1
          enddo
        enddo
      enddo
    enddo

  end subroutine update_eta

  subroutine update_vp_tiso()
    use global_var , only : FOUR_THIRDS, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220
    ! model update:
    ! transverse isotropic update only in layer Moho to 220
    ! (where SPECFEM3D_GLOBE considers TISO)
    ! everywhere else uses an isotropic update
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
      model_dbulk, model_dbetah, model_dbetav

    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL) :: alphav0, alphah0, betav0, betah0
    real(kind=CUSTOM_REAL) :: alphav1, alphah1
    real(kind=CUSTOM_REAL) :: dbetaiso, dbulk

    model_dbulk  = dmodels(:, :, :, :,KVPV_IDX)
    model_dbetav = dmodels(:, :, :, :, KVSV_IDX)
    model_dbetah = dmodels(:, :, :, :, KVSH_IDX)

    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! initial model values
            alphav0 = models(i,j,k,ispec,VPV_IDX)
            alphah0 = models(i,j,k,ispec,VPH_IDX)
            betav0  = models(i,j,k,ispec,VSV_IDX)
            betah0  = models(i,j,k,ispec,VSH_IDX)
            alphav1 = 0._CUSTOM_REAL
            alphah1 = 0._CUSTOM_REAL

            ! do not use transverse isotropy except if element is between d220 and Moho
            if(.not. ( idoubling(ispec)==IFLAG_670_220 .or. idoubling(ispec)==IFLAG_220_80 &
                       .or. idoubling(ispec)==IFLAG_80_MOHO) ) then

              ! isotropic model update
              ! shear values
              ! isotropic kernel K_beta = K_betav + K_betah
              ! with same scaling step_length the model update dbeta_iso = dbetav + dbetah
              ! note:
              !   this step length can be twice as big as that given by the input
              dbetaiso = model_dbetav(i,j,k,ispec) + model_dbetah(i,j,k,ispec)

              ! alpha values
              dbulk = model_dbulk(i,j,k,ispec)
              alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betav0**2 * ( &
                                  exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )

              alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betah0**2 * ( &
                                  exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
              ! note: alphah probably not used in SPECFEM3D_GLOBE
            else
              ! transverse isotropic model update
              ! alpha values
              dbulk = model_dbulk(i,j,k,ispec)
              alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) &
                              + FOUR_THIRDS * betav0**2 * ( &
                                  exp(2.0*model_dbetav(i,j,k,ispec)) - exp(2.0*dbulk) ) )

              alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) &
                              + FOUR_THIRDS * betah0**2 * ( &
                                  exp(2.0*model_dbetah(i,j,k,ispec)) - exp(2.0*dbulk) ) )

            endif
            ! stores new model values
            models_new(i,j,k,ispec,VPV_IDX) = alphav1
            models_new(i,j,k,ispec,VPH_IDX) = alphah1
          enddo
        enddo
      enddo
    enddo

  end subroutine update_vp_tiso

  subroutine update_vs_tiso()
    use global_var , only : FOUR_THIRDS, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220
    ! model update:
    ! transverse isotropic update only in layer Moho to 220
    ! (where SPECFEM3D_GLOBE considers TISO)
    ! everywhere else uses an isotropic update
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
      model_dbulk, model_dbetah, model_dbetav, model_deta

    integer :: i, j, k, ispec
    real(kind=CUSTOM_REAL) :: betav0, betah0
    real(kind=CUSTOM_REAL) :: betav1, betah1
    real(kind=CUSTOM_REAL) :: dbetaiso

    model_dbetav = dmodels(:, :, :, :, KVSV_IDX)
    model_dbetah = dmodels(:, :, :, :, KVSH_IDX)

    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

             ! initial model values
            betav0  = models(i,j,k,ispec,VSV_IDX)
            betah0  = models(i,j,k,ispec,VSH_IDX)
            betav1  = 0._CUSTOM_REAL
            betah1  = 0._CUSTOM_REAL


            ! do not use transverse isotropy except if element is between d220 and Moho
            if(.not. ( idoubling(ispec)==IFLAG_670_220 .or. idoubling(ispec)==IFLAG_220_80 &
                       .or. idoubling(ispec)==IFLAG_80_MOHO) ) then

              ! isotropic model update

              ! shear values
              ! isotropic kernel K_beta = K_betav + K_betah
              ! with same scaling step_length the model update dbeta_iso = dbetav + dbetah
              ! note:
              !   this step length can be twice as big as that given by the input
              dbetaiso = model_dbetav(i,j,k,ispec) + model_dbetah(i,j,k,ispec)
              betav1 = betav0 * exp( dbetaiso )
              betah1 = betah0 * exp( dbetaiso )

            else
              ! transverse isotropic model update
              ! shear values
              betav1 = betav0 * exp( model_dbetav(i,j,k,ispec) )
              betah1 = betah0 * exp( model_dbetah(i,j,k,ispec) )

            endif
            ! stores new model values
            models_new(i,j,k,ispec,VSV_IDX) = betav1
            models_new(i,j,k,ispec,VSH_IDX) = betah1

          enddo
        enddo
      enddo
    enddo

    
  end subroutine update_vs_tiso
  
 subroutine update_qmu()
    use global_var , only : FOUR_THIRDS, IFLAG_80_MOHO, IFLAG_220_80, IFLAG_670_220,QMU_IDX,KQMU_IDX
    
    ! model update:
    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: model_dqmu
    real(kind=CUSTOM_REAL) :: qmu0,qmu_min,qmu_max
    integer :: i, j, k, ispec


    qmu_min = 30.0
    qmu_max = 9000.0
    model_dqmu = dmodels(:,:,:,:,KQMU_IDX)

    do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX


              
              qmu0 = models(i,j,k,ispec,QMU_IDX)

           
                 
              models_new(i,j,k,ispec,QMU_IDX) = 1.0 / ( (1.0 / qmu0) * exp(model_dqmu(i,j,k,ispec)) )

              ! Limits for attenuation model
              if (models_new(i,j,k,ispec,QMU_IDX) > qmu_max) then
                 models_new(i,j,k,ispec,QMU_IDX) = qmu_max
              else if(models_new(i,j,k,ispec,QMU_IDX) < qmu_min) then
                 models_new(i,j,k,ispec,QMU_IDX) = qmu_min
              endif

           
             
           enddo
        enddo
     enddo
  enddo

  end subroutine update_qmu
  

  subroutine store_perturbations(outputfile)
    character(len=*),intent(in) :: outputfile

    real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB) :: models_perturb

    where(models /= 0.0)
      models_perturb = (models_new - models) / models
    elsewhere
      models_perturb = 0.0
    end where

    if(myrank == 0) print*, "Relative Model Update: "
    call stats_value_range(models_perturb, MODEL_PERTURB_NAMES_GLOB)

    call write_bp_file(models_perturb, MODEL_PERTURB_NAMES_GLOB, "PERTURBS_GROUP", outputfile)
  end subroutine store_perturbations

  subroutine save_output_files(outputdir)
    character(len=*), intent(in) :: outputdir

    character(len=500) :: outputfile

    ! stores new model in files
    outputfile = trim(outputdir)//'/model_gll.bp'
    if(myrank == 0) print*, "New model file: ", trim(outputfile)
    call write_bp_file(models_new, MODEL_NAMES_TOTAL , "KERNELS_GROUP", outputfile)

    outputfile = trim(outputdir)//'/dkernels.bp'
    if(myrank == 0) print*, "Kernel Change file: ", trim(outputfile)
    call write_bp_file(dmodels, KERNEL_NAMES_GLOB, "KERNELS_GROUP", outputfile)

    ! stores relative model perturbations
    if(myrank == 0) print*, "Model Perturbation file: ", trim(outputfile)
    outputfile = trim(outputdir)//'/perturbs_gll.bp'
    call store_perturbations(outputfile)
  end subroutine save_output_files

  subroutine stats_value_range(values, varnames)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: values
    character(len=*), dimension(:), intent(in) :: varnames

    integer :: n, i
    real(kind=CUSTOM_REAL) :: vmax, vmin

    n = size(varnames)

    do i=1,n
      call max_all_all_cr(maxval(values(:, :, :, :, i)), vmax)
      call min_all_all_cr(minval(values(:, :, :, :, i)), vmin)
      if(myrank == 0) then
        write(*, '(A30, A, ES16.8, A, ES16.8, A)') trim(varnames(i)), &
          " range -- ( ", vmin, ",", vmax, "  )"
      end if
    enddo

  end subroutine stats_value_range

end module model_update_par

! Program main
program main
  use mpi
  use adios_read_mod
  use model_update_par
  use global_var, only : myrank

  implicit none
  integer :: ier

  ! model update length
  real(kind=CUSTOM_REAL) :: step_fac
  ! program arguments, path of input and output files
  character(len=500) :: input_model_file, input_kernel_file, input_solver_file, outputdir,kernel_parfile

  call init_mpi()

  call get_sys_args(kernel_parfile,step_fac, input_model_file, input_solver_file, &
       input_kernel_file, outputdir)

  call init_module_par(kernel_parfile)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)
  ! reads in current transverse isotropic model files: vpv.. & vsv.. & eta & rho
  if(myrank == 0) print*, "|<---- Reading Model File ---->|"
  call read_bp_file_real(input_model_file, MODEL_NAMES_GLOB, models)

  if (NMODEL0 .gt. 0) then
  call read_bp_file_real(input_model_file, MODEL0_NAMES, models_fixed)
  endif

  ! reads in smoothed kernels: bulk, betav, betah, eta
  if(myrank == 0) print*, "|<---- Reading Kernel file ---->|"
  call read_bp_file_real(input_kernel_file, KERNEL_NAMES_GLOB, kernels)

  if(myrank == 0) print*, "|<---- Reading Solver(Mesh) File ---->|"
  call read_solver_data(input_solver_file)

  ! calculates gradient
  ! steepest descent method
  if(myrank == 0) print*, "|<---- Calculate Model Update ---->|"
  call get_model_change(step_fac)

  ! compute new model in terms of alpha, beta, eta and rho
  ! (see also Carl's Latex notes)
  if(myrank == 0) print*, "|<---- Apply Model Update --->|"
  call update_model_tiso()

  if(myrank == 0) print*, "|<---- Save Output --->|"
  call save_output_files(outputdir)

  if(myrank == 0) print*, "|<---- Done Writing ---->|"
  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program main
