module merge_files_mod
  use mpi
  use global_var, only : myrank, nprocs, NGLLX, NGLLY, NGLLZ, NSPEC, NGLOB, CUSTOM_REAL,exit_mpi
  use global_var, only : max_all_all_cr,min_all_all_cr
  implicit none


contains
  
    subroutine get_sys_args(config_file,output_file)
    ! input args
    character(len=*), intent(inout) :: config_file,output_file

    call getarg(1, config_file)
    call getarg(2, output_file)

    
    if (trim(config_file) == '' .or. trim(output_file) == '' ) then
        call exit_MPI('Usage: merge_adios_files config_file output_file')
    endif


    if (myrank == 0) then
      print*,'Merge adios files'
      print*, "System Args: "
      print*,' config_file  :', trim(config_file)
      print*,' output_file :', trim(output_file)
      print*
   endif
   
 end subroutine get_sys_args


 subroutine read_config_file(fn,file1,file1_par,npar1,file2,file2_par,npar2)
   character(len=*), intent(in) :: fn
   character(len=512), dimension(:), allocatable, intent(inout) :: file1_par, &
        file2_par
   character(len=*), intent(inout) :: file1,file2
   integer, intent(inout) :: npar1,npar2
   integer:: i
   integer :: fh = 1001

   if(myrank == 0) print*, '|<============= Parsing Config File =============>|'

    open(fh, file=trim(fn), status='old')

    read(fh, '(A)') file1

    if(myrank == 0) print*, "File 1: ",trim(file1)
    read(fh, *) npar1
    
    if(myrank == 0) print*, "number parameters model 1: ", npar1

    allocate(file1_par(npar1))

    do i=1,npar1
       read(fh, '(A)') file1_par(i)
       if(myrank == 0) write(*, '(A,A)') "Model file par 1: ", trim(file1_par(i))
    enddo


    read(fh, '(A)') file2

    if(myrank == 0) print*, "File 2: ", file2
    
    read(fh,*) npar2

    if(myrank == 0) print*, "number parameters model 2: ", npar2

    allocate(file2_par(npar2))
    
    do i=1,npar2
       read(fh, '(A)') file2_par(i)
       if(myrank == 0)  write(*, '(A,A)') "Model file par 2: " , trim(file2_par(i))
    enddo

    close(fh)
    
  end subroutine read_config_file


  subroutine merge_models(model1,model2,file1_par,file2_par,model,model_name)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: model1,model2
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: model
    character(len=*), dimension(:), intent(in) :: file1_par,file2_par
    character(len=512), dimension(:), allocatable, intent(inout) :: model_name
    integer :: npar_model,i

    npar_model = size(file1_par) + size(file2_par)

    if(myrank == 0) print*, "No. model parameters: ",npar_model
    allocate(model_name(npar_model))


    if(myrank == 0) print*, "No. model parameters <FILE 1> : ",size(file1_par)
    do i=1,size(file1_par)
       model(:,:,:,:,i) = model1(:,:,:,:,i)
       model_name(i) = file1_par(i)
    enddo

    if(myrank == 0) print*, "No. model parameters <FILE 2> : ",size(file2_par)
    do i=1,size(file2_par)
       model(:,:,:,:,i + size(file1_par)) = model2(:,:,:,:,i)
       model_name(i + size(file1_par)) = file2_par(i)
    enddo    
    
  end subroutine merge_models


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


  end module merge_files_mod

  ! Program main
program main
  use mpi
  use adios_read_mod
  use merge_files_mod
  use AdiosIO
  use global_var, only : myrank,NGLLX,NGLLY,NGLLZ,NSPEC,init_mpi

  implicit none
  integer :: ier,npar1,npar2,i
  integer :: fh=1001
  ! program arguments, path of input and output files
  character(len=500) :: file1,file2,config_file,output_file
  character(len=500), dimension(:),allocatable :: file1_par,file2_par,model_name
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: model1,model2,model

  call init_mpi()

  call get_sys_args(config_file,output_file)

  if(myrank == 0) print*, "|<---- Reading Config File ---->|"
!  call read_config_file(config_file,file1,file1_par,npar1,file2,file2_par,npar2)

  if(myrank == 0) print*, '|<============= Parsing Config File =============>|'

    open(fh, file=trim(config_file), status='old')
    read(fh, '(A)') file1
    if(myrank == 0) print*, "File 1: ",trim(file1)
    read(fh, *) npar1
    if(myrank == 0) print*, "number parameters model 1: ", npar1

    allocate(file1_par(npar1))

    do i=1,npar1
       read(fh, '(A)') file1_par(i)
       if(myrank == 0) write(*, '(A,A)') "Model file par 1: ", trim(file1_par(i))
    enddo


    read(fh, '(A)') file2
    if(myrank == 0) print*, "File 2: ", trim(file2)
    read(fh,*) npar2
    if(myrank == 0) print*, "number parameters model 2: ", npar2

    allocate(file2_par(npar2))
    
    do i=1,npar2
       read(fh, '(A)') file2_par(i)
       if(myrank == 0)  write(*, '(A,A)') "Model file par 2: " , trim(file2_par(i))
    enddo

    close(fh)

  if(myrank == 0) print*, "|<---- Allocating memory for models ---->|"

  allocate(model1(NGLLX,NGLLY,NGLLZ,NSPEC,npar1),model2(NGLLX,NGLLY,NGLLZ,NSPEC,npar2),&
       model(NGLLX,NGLLY,NGLLZ,NSPEC,npar1 + npar2))
  
  call adios_read_init_method (ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
       "verbose=1", ier)
  
  ! reads in current transverse isotropic model files: vpv.. & vsv.. & eta & rho
  if(myrank == 0) print*, "|<---- Reading Model File 1 ---->|"
  call read_bp_file_real(file1, file1_par, model1)
  call stats_value_range(model1,file1_par)

  if(myrank == 0) print*, "|<---- Reading Model File 2 ---->|"
  call read_bp_file_real(file2, file2_par, model2)
  call stats_value_range(model2,file2_par)

  if(myrank == 0) print*, "|<---- Merging Model Files  ---->|"
  call merge_models(model1,model2,file1_par,file2_par,model,model_name)


  call stats_value_range(model,model_name)

  
  if(myrank == 0) print*, "|<---- Save Output --->|"
  call write_bp_file(model,model_name,"KERNELS_GROUP",output_file)

  if(myrank == 0) print*, "|<---- Done Writing ---->|"
  call adios_finalize(myrank, ier)
  call MPI_FINALIZE(ier)

end program main
