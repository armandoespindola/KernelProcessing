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
  
  
contains


  subroutine get_args(kernel_parfile,solver_file,input_model_file,perturb_file,output_dir)
    character(len=*), intent(inout) :: output_dir,solver_file,input_model_file
    character(len=*), intent(inout) :: kernel_parfile,perturb_file

    call getarg(1, kernel_parfile)
    call getarg(2, solver_file)
    call getarg(3, input_model_file)
    call getarg(4, perturb_file)
    call getarg(5, output_dir)

    if(trim(output_dir) == '' .or. trim(kernel_parfile) == '' .or. trim(input_model_file) == '' &
         .or. trim(perturb_file) == '' .or. trim(solver_file) == '') then
      call exit_mpi("Usage: xperturb_model kernel_parfile solver_file input_model_file perturb_file output_dir")
    endif

    if(myrank == 0) then
       write(*, *) "kernel parfile: ", trim(kernel_parfile)
       write (*,*) "perturb file: ", trim(perturb_file)
       write (*,*) "input model file: ", trim(input_model_file)
       write(*, *) "Solver file: ", trim(solver_file)
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

  character(len=500) :: input_model_file,output_dir,solver_file,output_file,perturb_file
  character(len=500) :: kernel_parfile,kernel_name,model_name
  integer :: ier,i,idx_per,np,idx_model,j,k,ispec

  real(kind=CUSTOM_REAL),dimension(:),allocatable :: depth,lat,sigmav,sigmah,lon,fact
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: models
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: models_new,perturbs

  integer :: fh = 1001
  character(len=500) :: line


  call init_mpi()


  call get_args(kernel_parfile,solver_file,input_model_file,perturb_file,output_dir)
  call init_kernel_par(kernel_parfile)

  allocate(models(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB),models_new(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB))
  allocate(perturbs(NGLLX,NGLLY,NGLLZ,NSPEC,NKERNEL_GLOB))
  
  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
       "verbose=1", ier)



  if(myrank == 0) write(*, '(A)') "|<-------- Reading solver file -------->|"
  call read_model_file(solver_file)

  if(myrank == 0) write(*, '(A)') "|<-------- Reading file perturbations -------->|"
  open(fh,file=trim(perturb_file), status='old')
  read(fh, '(A)') line
  kernel_name = trim(line)
  if(myrank == 0) print*, "Kernel Name: ",trim(kernel_name)
  read(fh, '(A)') line
  model_name = trim(line)
  if(myrank == 0) print*, "Model Name: ",trim(model_name)
  read(fh,*) np

  allocate(lat(np),lon(np),depth(np),sigmav(np),sigmah(np),fact(np))
  do i=1,np
     read(fh, *) lat(i),lon(i),depth(i),sigmav(i),sigmah(i),fact(i)
  enddo

  close(fh)


  do i=1,NKERNEL_GLOB
     if (trim(kernel_name) == trim(KERNEL_NAMES_GLOB(i))) then
        idx_per = i
     endif
  enddo

  do i=1,NPAR_GLOB
     if ("reg1/"//trim(model_name) == trim(MODEL_NAMES_GLOB(i))) then
        idx_model = i
     endif
  enddo


  if(myrank == 0) print*, "|<---- Reading Model File ---->|"
  call read_bp_file_real(input_model_file, MODEL_NAMES_GLOB, models)
  
  perturbs(:,:,:,:, :) = 0.0_CUSTOM_REAL

  if(myrank == 0) write(*, '(A)') "|<-------- Creating perturbations -------->|"

  do i=1,np
     if (myrank==0) write (*,*) "r, lat, lon, sigmah0, sigmav0, fact"
     if (myrank==0) write (*,*)  depth(i), lat(i),lon(i), sigmah(i), sigmav(i),fact(i)

     call add_gaussian_perturb_elliptic(depth(i), lat(i), &
          lon(i), sigmah(i), sigmav(i),idx_per,perturbs)
     
     perturbs(:,:,:,:,idx_per) = perturbs(:,:,:,:,idx_per) +  perturbs(:,:,:,:,idx_per) * fact(i)
     
  enddo


  models_new = models
  
  if(myrank == 0) write(*, '(A)') "|<-------- Adding perturbations -------->|"


  
  if (trim(model_name) == "qmu") then
     do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              models_new(i,j,k,ispec,idx_model) = 1.0 / ( (1.0 / models_new(i,j,k,ispec,idx_model)) * &
                   exp(perturbs(i,j,k,ispec,idx_per)))
           enddo
        enddo
     enddo
  enddo
  
else
   do ispec = 1, NSPEC
      do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              models_new(i,j,k,ispec,idx_model) = models_new(i,j,k,ispec,idx_model) * exp(perturbs(i,j,k,ispec,idx_per))
           enddo
        enddo
     enddo
  enddo
  endif
  
  
     
  ! stores new model in files
  output_file = trim(output_dir)//'/perturb_m.bp'
  if(myrank == 0) print*, "Perturb file: ", trim(output_file)
  call write_bp_file(perturbs, KERNEL_NAMES_GLOB, "KERNELS_GROUP", output_file)


  output_file = trim(output_dir)//'/model_gll_dm.bp'
  if(myrank == 0) print*, "Model file: ", trim(output_file)
  call write_bp_file(models_new, MODEL_NAMES_GLOB, "KERNELS_GROUP", output_file)

 call adios_finalize(myrank, ier)
 call MPI_FINALIZE(ier)
  

  
  
  end program 
