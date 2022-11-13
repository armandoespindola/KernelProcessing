module global_var

  use mpi
  implicit none

  integer, parameter :: CUSTOM_REAL = 4
  integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL

  ! values from the mesh
  double precision, parameter :: GAUSSALPHA = 0.d0, GAUSSBETA = 0.d0
  integer, parameter :: NGLLX = 5, NGLLY = NGLLX, NGLLZ = NGLLX
  integer, parameter :: NSPEC_CRUST_MANTLE = 492
  integer, parameter :: NGLOB_CRUST_MANTLE = 35573
  integer, parameter :: NSPEC = NSPEC_CRUST_MANTLE
  integer, parameter :: NGLOB = NGLOB_CRUST_MANTLE

  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD = "MPI"
  integer, parameter :: ADIOS_BUFFER_SIZE_IN_MB = 1000

  ! Earth
  double precision, parameter :: R_EARTH = 6371000.d0
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0

  real(kind=CUSTOM_REAL), parameter :: FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: TWO_THIRDS = 2._CUSTOM_REAL/3._CUSTOM_REAL
  real(kind=CUSTOM_REAL), parameter :: ONE_THIRD = 1._CUSTOM_REAL/3._CUSTOM_REAL

  ! flags for idoubling
  integer, parameter :: IFLAG_80_MOHO = 2
  integer, parameter :: IFLAG_220_80 = 3
  integer, parameter :: IFLAG_670_220 = 4

  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0

  ! MPI variable that need to be init in the main program
  integer :: myrank, nprocs

  !parameters for kernel processing
  integer :: NPAR_GLOB,NKERNEL_GLOB,NHESS0
  integer :: QMU_IDX,KQMU_IDX,NMODEL0,NMODEL_TOTAL
  logical :: ATTENUATION_FLAG
  character(len=500), dimension(:),allocatable :: KERNEL_NAMES_GLOB,MODEL_NAMES_GLOB,MODEL_PERTURB_NAMES_GLOB
  character(len=500), dimension(:),allocatable :: KER_HESS_NAMES_GLOB,HESS_NAMES_GLOB
  character(len=500), dimension(:),allocatable :: MODEL0_NAMES,MODEL_NAMES_TOTAL,KERNEL_NAMES_GLOB_NQ
  character(len=500) :: MODEL
  integer :: VPV_IDX,VPH_IDX,VSV_IDX,VSH_IDX,RHO_IDX,ETA_IDX,VP_IDX,VS_IDX
  integer :: KVPV_IDX,KVPH_IDX,KVSV_IDX,KVSH_IDX,KRHO_IDX,KETA_IDX,KVP_IDX,KVS_IDX
  

  contains

  subroutine sum_all_all_cr(sendbuf, recvbuf)
    real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
    integer :: ier

    call MPI_AllReduce(sendbuf, recvbuf, 1, CUSTOM_MPI_TYPE, MPI_SUM, MPI_COMM_WORLD, ier)
  end subroutine

  subroutine max_all_all_cr(sendbuf, recvbuf)

    real(kind=CUSTOM_REAL):: sendbuf, recvbuf
    integer :: ier

    call MPI_ALLREDUCE(sendbuf, recvbuf, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD, ier)

  end subroutine max_all_all_cr

  subroutine min_all_all_cr(sendbuf, recvbuf)

    real(kind=CUSTOM_REAL):: sendbuf, recvbuf
    integer :: ier

    call MPI_ALLREDUCE(sendbuf, recvbuf, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD, ier)

  end subroutine min_all_all_cr

  subroutine build_jacobian(ibool, xix, xiy, xiz, etax, etay, etaz, &
                            gammax, gammay, gammaz, jacobian)
  integer, dimension(:, :, :, :), intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: xix, xiy, xiz
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: etax, etay, etaz
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: gammax, gammay, gammaz
  real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(inout) :: jacobian

  integer :: i, j, k, ispec, iglob
  real(kind=CUSTOM_REAL) :: xixl, xiyl, xizl, etaxl, etayl, etazl
  real(kind=CUSTOM_REAL) :: gammaxl, gammayl, gammazl, jacobianl

  do ispec = 1, NSPEC
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! build jacobian
          ! get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          ! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))
          jacobian(i, j, k, ispec) = jacobianl
        enddo
      enddo
    enddo
  enddo

  end subroutine build_jacobian

  subroutine build_gll_weight(wgll_cube)
    double precision, dimension(:, :, :), intent(inout) :: wgll_cube

    ! Gauss-Lobatto-Legendre points of integration and weights
    double precision, dimension(NGLLX) :: xigll, wxgll
    double precision, dimension(NGLLY) :: yigll, wygll
    double precision, dimension(NGLLZ) :: zigll, wzgll

    integer :: i, j, k

    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
           wgll_cube(i,j,k) = wxgll(i)*wygll(j)*wzgll(k)
        enddo
      enddo
    enddo
  end subroutine build_gll_weight

  subroutine Parallel_ComputeInnerProduct(vect1, vect2, Niv, jacobian, qp)
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in) :: vect1, vect2
    integer,                intent(in)    :: Niv
    real(kind=CUSTOM_REAL), dimension(:, :, :, :) :: jacobian
    real(kind=CUSTOM_REAL), intent(inout) :: qp

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, Niv) :: wks_1n, wks_2n

    real(kind=CUSTOM_REAL)  :: qp_tmp_single
    ! try double precision
    real(kind=8)            :: jacobianl, weight, qp_tmp
    integer                 :: ipar, i, j, k, ispec
    real(kind=CUSTOM_REAL)  :: coeff, coeff_n1, coeff_n2
    real(kind=8)            :: coeff_n1_dp, coeff_n2_dp

    double precision, dimension(NGLLX, NGLLY, NGLLZ) :: wgll_cube

    !! try normalization to avoid numerical errors
    !call Parallel_ComputeL2normSquare(vect1 , Niv, coeff_n1)
    !call Parallel_ComputeL2normSquare(vect2 , Niv, coeff_n2)

    coeff=maxval(abs(vect1(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n1)
    if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL
    wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1

    coeff=maxval(abs(vect2(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n2)
    if (coeff_n2 == 0._CUSTOM_REAL) coeff_n2=1._CUSTOM_REAL
    wks_2n(:,:,:,:,:) = vect2(:,:,:,:,:) / coeff_n2

    coeff_n1_dp = coeff_n1
    coeff_n2_dp = coeff_n2

    call build_gll_weight(wgll_cube)

    qp_tmp=0._CUSTOM_REAL
    do ipar=1, Niv
      do ispec = 1, NSPEC
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              weight = wgll_cube(i, j, k)
              jacobianl = jacobian(i, j, k, ispec)
              qp_tmp = qp_tmp + jacobianl * weight * &
                       wks_1n(i,j,k,ispec,ipar) * wks_2n(i,j,k,ispec,ipar)
              !qp = qp + jacobianl * weight * &
              !          vect1(i,j,k,ispec,ipar) * vect2(i,j,k,ispec,ipar)
            enddo
          enddo
        enddo
      enddo
    enddo

    qp_tmp_single = REAL(qp_tmp * coeff_n1_dp * coeff_n2_dp, CUSTOM_REAL)
    qp=0.
    call sum_all_all_cr(qp_tmp_single, qp)

  end subroutine Parallel_ComputeInnerProduct

  subroutine Parallel_ComputeL2normSquare(vect1 , Niv, jacobian, qp)
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), intent(in) :: vect1
    integer,                intent(in)    :: Niv
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: jacobian
    real(kind=CUSTOM_REAL), intent(inout) :: qp

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, Niv) :: wks_1n
    real(kind=CUSTOM_REAL) :: coeff, coeff_n1
    real(kind=CUSTOM_REAL) :: qp_tmp
    real(kind=8)           :: jacobianl, weight, qp_dp, coeff_n1_dp
    integer                :: ipar, i, j, k, ispec

    ! Gauss-Lobatto-Legendre points of integration and weights
    double precision, dimension(NGLLX, NGLLY, NGLLZ) :: wgll_cube

    coeff=maxval(abs(vect1(:,:,:,:,:)))
    call max_all_all_cr(coeff, coeff_n1)

    if (coeff_n1 == 0._CUSTOM_REAL) coeff_n1=1._CUSTOM_REAL

    wks_1n(:,:,:,:,:) = vect1(:,:,:,:,:) / coeff_n1
    coeff_n1_dp=coeff_n1

    call build_gll_weight(wgll_cube)

    qp_dp=0.d0
    do ipar=1,Niv
      do ispec = 1, NSPEC
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              weight = wgll_cube(i, j, k)
              jacobianl = jacobian(i, j, k, ispec)
              qp_dp = qp_dp + jacobianl * weight * (wks_1n(i,j,k,ispec,ipar)**2)
            enddo
           enddo
        enddo
      enddo
    enddo

    qp_tmp = REAL(qp_dp * coeff_n1_dp * coeff_n1_dp, CUSTOM_REAL)
    qp = 0.d0
    call sum_all_all_cr(qp_tmp, qp)

  end subroutine Parallel_ComputeL2normSquare

  subroutine Parallel_ComputeIntegral(vec, jacobian, res)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: vec
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian
    real(kind=CUSTOM_REAL), intent(inout) :: res

    double precision, dimension(NGLLX, NGLLY, NGLLZ) :: wgll_cube
    integer :: i, j, k, ispec
    real(kind=8) :: integl, weight, jacobianl
    real(kind=4) :: integl_tmp

    call build_gll_weight(wgll_cube)

    integl = 0.0
    do ispec = 1, NSPEC
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            weight = wgll_cube(i, j, k)
            jacobianl = jacobian(i, j, k, ispec)
            integl = integl + jacobianl * weight * vec(i,j,k,ispec)
          enddo
         enddo
      enddo
    enddo

    integl_tmp = REAL(integl, CUSTOM_REAL)
    res = 0.0
    call sum_all_all_cr(integl_tmp, res)

  end subroutine Parallel_ComputeIntegral

  subroutine exit_mpi(error_msg)
    character(len=*) :: error_msg

    integer :: ier

    
    ! write error message to screen
    write(*,*) error_msg(1:len(error_msg))
    write(*,*) 'Error detected, aborting MPI... proc ',myrank
    call MPI_ABORT(MPI_COMM_WORLD, -1, ier)
  end subroutine exit_mpi

  subroutine init_mpi()
    integer :: ier
    call MPI_INIT(ier)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ier)
    if (myrank == 0) print*, "MPI init with nprocs=", nprocs
  end subroutine init_mpi

  subroutine init_kernel_par(kernel_parfile)
    character(len=*), intent(inout) :: kernel_parfile

    integer :: fh = 1001
    integer :: i,ier,k
    character(len=500) :: line
    ATTENUATION_FLAG = .false.

    VPV_IDX = -1
    VPH_IDX = -1
    VSV_IDX = -1
    VSH_IDX = -1
    RHO_IDX = -1
    ETA_IDX = -1
    VP_IDX = -1
    VS_IDX = -1

    
    KVPV_IDX = -1
    KVPH_IDX = -1
    KVSV_IDX = -1
    KVSH_IDX = -1
    KRHO_IDX = -1
    KETA_IDX = -1
    KVP_IDX = -1
    KVS_IDX = -1
    
    ! NPAR_GLOB=7
    ! NKERNEL_GLOB=1

    ! allocate(KERNEL_NAMES_GLOB(NKERNEL_GLOB),MODEL_NAMES_GLOB(NPAR_GLOB),MODEL_PERTURB_NAMES_GLOB(NPAR_GLOB),stat=ier)
    
    
    ! MODEL_NAMES_GLOB(1)="reg1/rho"
    ! MODEL_NAMES_GLOB(2)="reg1/vpv"
    ! MODEL_NAMES_GLOB(3)="reg1/vph"
    ! MODEL_NAMES_GLOB(4)="reg1/vsv"
    ! MODEL_NAMES_GLOB(5)="reg1/vsh"
    ! MODEL_NAMES_GLOB(6)="reg1/eta"
    ! MODEL_NAMES_GLOB(7)="reg1/qmu"


    ! MODEL_PERTURB_NAMES_GLOB(1)="reg1/drhorho"
    ! MODEL_PERTURB_NAMES_GLOB(2)="reg1/dvpvvpv"
    ! MODEL_PERTURB_NAMES_GLOB(3)="reg1/dvphvph"
    ! MODEL_PERTURB_NAMES_GLOB(4)="reg1/dvsvvsv"
    ! MODEL_PERTURB_NAMES_GLOB(5)="reg1/dvshvsh"
    ! MODEL_PERTURB_NAMES_GLOB(6)="reg1/detaeta"
    ! MODEL_PERTURB_NAMES_GLOB(7)="reg1/dqmuqmu"

    ! KERNEL_NAMES_GLOB(1)="mu_kl_crust_mantle"

    ! QMU_IDX = 7
    ! KQMU_IDX = 1
    ! ATTENUATION_FLAG=.true.
    

    !OPENING FILE TO READ PARAMETERS

    if (myrank == 0) print*, "KERNEL PARFILE: ",trim(kernel_parfile)

    open(fh,file=trim(kernel_parfile), status='old')
    read(fh, '(A)') line
    MODEL = trim(line)
    read(fh,*) NPAR_GLOB
    !if (myrank == 0) print*, NPAR_GLOB
    
    allocate(MODEL_NAMES_GLOB(NPAR_GLOB), &
         MODEL_PERTURB_NAMES_GLOB(NPAR_GLOB),stat=ier)

    do i=1,NPAR_GLOB
       read(fh, '(A)') line
       !if (myrank == 0) print*,trim(line)
       MODEL_NAMES_GLOB(i)="reg1/"//trim(line)
       MODEL_PERTURB_NAMES_GLOB(i)="reg1/d"//trim(line)//"d"//trim(line)
    end do

    read(fh,*) NKERNEL_GLOB
    !if (myrank == 0) print*,NKERNEL_GLOB

    allocate(KERNEL_NAMES_GLOB(NKERNEL_GLOB))
    
    do i=1,NKERNEL_GLOB
       read(fh, '(A)') line
       !if (myrank == 0) print*,trim(line)
       KERNEL_NAMES_GLOB(i)=trim(line)
       KER_HESS_NAMES_GLOB(i)=trim(line)
    !   KER_HESS_NAMES_GLOB(i + NKERNEL_GLOB)="hess_"//trim(line)
    !   HESS_NAMES_GLOB(i)="hess_"//trim(line)
    end do

    read(fh,*) NMODEL0
    !if (myrank == 0) print*,NMODEL0

    if (NMODEL0 .gt. 0) then
       allocate(MODEL0_NAMES(NMODEL0),stat=ier)
       do i=1,NMODEL0
          read(fh, '(A)') line
          !if (myrank == 0) print*,trim(line)
          MODEL0_NAMES(i)="reg1/"//trim(line)
       end do
    endif

    read(fh,*) NHESS0

    allocate(KER_HESS_NAMES_GLOB(NKERNEL_GLOB + NHESS0),stat=ier)
    allocate(HESS_NAMES_GLOB(NHESS0),stat=ier)
    
    do i=1,NHESS0
       read(fh, '(A)') line
       !if (myrank == 0) print*,trim(line)
       KER_HESS_NAMES_GLOB(i + NKERNEL_GLOB)=trim(line)
       HESS_NAMES_GLOB(i)=trim(line)
    end do
    
    
    
    close(fh)

    NMODEL_TOTAL = NPAR_GLOB + NMODEL0
    allocate(MODEL_NAMES_TOTAL(NMODEL_TOTAL),stat=ier)

    do i=1,NPAR_GLOB
       MODEL_NAMES_TOTAL(i) = trim(MODEL_NAMES_GLOB(i))
    enddo

    if (NMODEL0 .gt. 0) then
       do i=NPAR_GLOB+1,NMODEL_TOTAL
          MODEL_NAMES_TOTAL(i) = trim(MODEL0_NAMES(i - NPAR_GLOB))
       enddo
    end if



    if (trim(MODEL) == "TISO") then
       VPV_IDX = check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/vpv")
       VPH_IDX = check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/vph")
       VSV_IDX = check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/vsv")
       VSH_IDX = check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/vsh")
       RHO_IDX = check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/rho")
       ETA_IDX = check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/eta")
       QMU_IDX= check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/qmu")

       KVPV_IDX = check_model(KERNEL_NAMES_GLOB,NKERNEL_GLOB,"bulk_c_kl_crust_mantle")
       KVPH_IDX = KVPV_IDX 
       KVSV_IDX = check_model(KERNEL_NAMES_GLOB,NKERNEL_GLOB,"bulk_betav_kl_crust_mantle")
       KVSH_IDX = check_model(KERNEL_NAMES_GLOB,NKERNEL_GLOB,"bulk_betah_kl_crust_mantle")
       KETA_IDX = check_model(KERNEL_NAMES_GLOB,NKERNEL_GLOB,"eta_kl_crust_mantle")
    endif

    if (trim(MODEL) == "ISO") then
       VP_IDX =  check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/vp")
       VS_IDX =  check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/vs")
       RHO_IDX =  check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/rho")
       QMU_IDX= check_model(MODEL_NAMES_GLOB,NPAR_GLOB,"reg1/qmu")
       
       KVP_IDX = check_model(KERNEL_NAMES_GLOB,NKERNEL_GLOB,"bulk_c_kl_crust_mantle")
       KVS_IDX = check_model(KERNEL_NAMES_GLOB,NKERNEL_GLOB,"bulk_beta_kl_crust_mantle")
    
    endif
    
 

    if (QMU_IDX .gt. 0 ) then
       ATTENUATION_FLAG=.true.
       KQMU_IDX= check_model(KERNEL_NAMES_GLOB,NKERNEL_GLOB,"mu_kl_crust_mantle")
    endif
    

    if ((ATTENUATION_FLAG) .and. (NPAR_GLOB .gt. 1)) then
       ! KERNEL_NAMES_GLOB WIHOUT ATTENUATION_NAME
       allocate(KERNEL_NAMES_GLOB_NQ(NKERNEL_GLOB - 1))
       k=1
       do i=1,NKERNEL_GLOB
          if ("mu_kl_crust_mantle" .ne. trim(KERNEL_NAMES_GLOB(i))) then
             KERNEL_NAMES_GLOB_NQ(k) = trim(KERNEL_NAMES_GLOB(i))
             k = k + 1
          endif
          
       enddo
       
    endif


    ! if (myrank ==  0) print*, "########## Number of parameters for Kernel Processing #########"
    ! if (myrank ==  0) print*, "Npar: ",NPAR_GLOB

    ! allocate(KERNEL_NAMES_GLOB(NKERNEL_GLOB),MODEL_NAMES_GLOB(NPAR_GLOB),UPDATE_FLAG_GLOB(NPAR_GLOB), &
    !      SCALING_GLOB(NPAR_GLOB),MODEL_PERTURB_NAMES_GLOB(NKERNEL_GLOB),stat=ier)

    ! KERNEL_NAMES_GLOB=""
    ! MODEL_NAMES_GLOB="reg1/"
    ! MODEL_PERTURB_NAMES_GLOB="reg1/"
    ! UPDATE_FLAG_GLOB=0
    ! SCALING_GLOB=0.0d0

    ! if (ier /= 0) stop " Error Allocating parameter for Kernel Processing "



    ! do i=1,NPAR_GLOB
    !    read(fh, '(A)') line
    !   call get_par(line,i)
    ! end do


    
    if (myrank==0) then
       do i=1,NPAR_GLOB
          write(*, '(A,A)')"MODEL NAMES: " , trim(MODEL_NAMES_GLOB(i))
       enddo
    endif

   
    if (myrank==0) then
       do i=1,NKERNEL_GLOB
          write(*, '(A,A,X,A)')"KERNEL NAMES: " , trim(KERNEL_NAMES_GLOB(i)),trim(MODEL_PERTURB_NAMES_GLOB(i))
       enddo
    endif

    if (myrank==0) then
       do i=1,NMODEL0
          write(*, '(A,A)')"FIXED MODEL NAMES: " , trim(MODEL0_NAMES(i))
       enddo
    endif

    if (myrank==0) then
       do i=1,NHESS0
          write(*, '(A,A)')"FIXED MODEL NAMES: " , trim(HESS_NAMES_GLOB(i))
       enddo
    endif

    
    ! close(fh)
    
  end subroutine init_kernel_par


  function check_model(varlist,nvarlist,var_name) result(flag)
    character(len=*), dimension(:), intent(in) :: varlist
    character(len=*), intent(in) :: var_name
    integer , intent (in) :: nvarlist
    integer :: flag,i
    flag=-1
    do i=1,nvarlist
       if (trim(varlist(i)) == trim(var_name)) then
          flag = i
          if (myrank==0) write(*,'(A, I3)')trim(var_name),flag
          exit
       endif
    enddo  

  end function check_model
    

  ! subroutine get_par(line,iline)
  !   character(len=*),intent(inout) :: line
  !   integer,intent(inout) :: iline
  !   integer :: idx0,idx1,offset,par,len_str
  !   character(len=500),dimension(4) :: par_data=""
    

  !   idx0=1
  !   idx1=len(line)
  !   par=1

  !   do while (index(line(idx0:),",") > 0 )
  !      offset = index(line(idx0:),",")
  !      idx1 = offset + idx0
  !      !print*,line(idx0:idx1-2),idx0,idx1
  !      par_data(par)(:len(line(idx0:idx1-1))) = line(idx0:idx1-2)
  !      idx0 = idx1
  !      par= par + 1
  !   end do

  !   if (len(trim(line(idx0:))) > 0 ) then
  !      !print*, "last ", (line(idx0-1:))       
  !      par_data(par)(:len(trim(line(idx0:))))=trim(line(idx0:))
  !   endif

  !   ! model_name

  !   len_str=len(trim(par_data(1)))
  !   MODEL_NAMES_GLOB(iline)=trim(MODEL_NAMES_GLOB(iline)) // trim(par_data(1))
    
  !   MODEL_PERTURB_NAMES_GLOB(iline)=trim(MODEL_PERTURB_NAMES_GLOB(iline)) // "d" // trim(par_data(1)) // &
  !        trim(par_data(1))
    
  !   ! kernel_name
  !   len_str=len(trim(par_data(3)))
  !   KERNEL_NAMES_GLOB(iline)=trim(par_data(3))

  !   ! update model flag

  !   if (len(trim(par_data(2))) > 0) then
  !      read(par_data(2),*)UPDATE_FLAG_GLOB(iline)
  !   else
  !      UPDATE_FLAG_GLOB(iline)=0
  !   endif

  !   if (len(trim(par_data(4))) > 0) then
  !      read(par_data(4),*) SCALING_GLOB(iline)
  !   else
  !      SCALING_GLOB(iline)=0
  !   endif
    
  !   !

  !   ! scaling

  !   !read(par_data(4),*)scaling
    
    
    
  ! end subroutine get_par
  

end module global_var
