  ! Module with parameters

module perturb_model_parameters
  use mpi
  use adios_read_mod
  use global_var


  ! ======================================================
  ! MODELS
  !integer, parameter :: NMODELS = 4
  !character(len=500), dimension(NMODELS), parameter :: model_names = &
  !  (/character(len=500) :: "reg1/vp", "reg1/vs", "reg1/rho", &
  !                          "reg1/qmu" /)
  ! MODEL ARRAY
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB) :: models = 0.0
  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB) :: models_new = 0.0

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: models
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: models_new

  
  
contains


  subroutine get_args(kernel_parfile,input_file,input_solver_file, output_dir)
    character(len=*), intent(inout) :: input_file, input_solver_file,output_dir
    character(len=*), intent(inout) :: kernel_parfile

    call getarg(1, kernel_parfile)
    call getarg(2, input_file)
    call getarg(3,input_solver_file)
    call getarg(4, output_dir)

    if(trim(input_solver_file) == '' .or. trim(input_file) == ''  &
         .or. trim(output_dir) == '' .or. trim(kernel_parfile) == '') then
      call exit_mpi("Usage: xperturb_model input_file input_solverfile output_dir")
    endif

    if(myrank == 0) then
       write(*, *) "Input ADIOS model: ", trim(input_file)
       write (*,*) "Input ADIOS solver file: ", trim(input_solver_file) 
      write(*, *) "Output dir: ", trim(output_dir)
   endif
   
  end subroutine get_args

  

end module perturb_model_parameters

program main
  use mpi
  use perturb_model_parameters
  use adios_read_mod
  use global_var, only: myrank,CUSTOM_REAL,init_kernel_par,NGLLX,NGLLY,NGLLZ,NPAR_GLOB,MODEL_NAMES_GLOB,NSPEC
  use gpp_utils

  implicit none

  character(len=500) :: input_file,output_dir,input_solver_file,output_file
  character(len=500) :: kernel_parfile
  integer :: ier

  real(kind=CUSTOM_REAL),dimension(1) :: depth,lat,sigmav,sigmah
  real(kind=CUSTOM_REAL),dimension(6) :: lon
  real(kind=CUSTOM_REAL):: fact
  integer::iz,ilat,ilon,idx_par


  call init_mpi()


  call get_args(kernel_parfile,input_file,input_solver_file, output_dir)
  call init_kernel_par(kernel_parfile)

  allocate(models(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB),models_new(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB))
  
  call init_vars()
  
  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  if (myrank==0) print*, "### READING MODEL PARAMETERS ####"
  call read_bp_file_real(input_file,MODEL_NAMES_GLOB,models)


  if (myrank==0) print*, "### READING SOLVER FILE ####"
  call read_model_file(input_solver_file)


  if(myrank == 0) write(*, '(A)') "|<-------- Add gaussian perturb -------->|"
  kernels(:,:,:,:, :) = 0.0_CUSTOM_REAL

  if(myrank == 0) write(*, '(A)') "|<-------- adding perturbations -------->|"

 




  ! Perturbations only attenuation
  ! idx_par = 4 -> dqu
  idx_par = 4
  models_new = models
  models_new(:,:,:,:,idx_par) = 0.0d0
  !depth=(/200.0,400.0,600.0 /)
  depth=(/400.0/)
  lat=(/0.0/)
  lon=(/-135.0,-45.0,45.0,0.0,135.0,180.0/)
  sigmav=(/350.0/)
  sigmah=(/1250.0/)

  do iz=1,1
     do ilat=1,1
        do ilon=1,6
           !add_gaussian_perturb_elliptic(r, lat, lon, sigmah0, sigmav0, perturb_idx)
          if (myrank==0) write (*,*) "r, lat, lon, sigmah0, sigmav0, perturb_idx"
          if (myrank==0) write (*,*)  depth(iz), lat(ilat),lon(ilon), sigmah(iz), sigmav(iz)

          

         !  if (iz == 2) then
       !    call add_gaussian_perturb_elliptic(depth(iz), lat(ilat), &
       !         lon(ilon) + 22.0, sigmah(iz), sigmav(iz),1)
       ! else

          call add_gaussian_perturb_elliptic(depth(iz), lat(ilat), &
                lon(ilon), sigmah(iz), sigmav(iz),1)
          
          !endif                 
           
          fact = 1.0 !(-1.0)**(iz + ilon + ilat)
          !kernels=0.0d0
          models_new(:,:,:,:,idx_par) = models_new(:,:,:,:,idx_par) +  kernels(:,:,:,:,1) * fact
          kernels= 0.0d0
        enddo
     enddo
     enddo

     models_new(:,:,:,:,idx_par) = 1.0/ ( (1.0 / models(:,:,:,:,idx_par)) * (1.0 + 0.65 * models_new(:,:,:,:,idx_par)) )
    
     
  ! stores new model in files
  output_file = trim(output_dir)//'/model_gll.bp'
  if(myrank == 0) print*, "New model file: ", trim(output_file)
 call write_bp_file(models_new, MODEL_NAMES_GLOB, "KERNELS_GROUP", output_file)


 call adios_finalize(myrank, ier)
 call clean_vars()
  call MPI_FINALIZE(ier)
  

  
  
  end program 
