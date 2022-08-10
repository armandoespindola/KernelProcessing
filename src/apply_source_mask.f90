module src_mask_subs

  use global_var, only : CUSTOM_REAL, myrank, exit_mpi
  use global_var, only : KERNEL_NAMES_GLOB,MODEL_NAMES_GLOB,MODEL_PERTURB_NAMES_GLOB
  use global_var, only : init_kernel_par,NPAR_GLOB,NKERNEL_GLOB
  implicit none

  contains

  subroutine create_name_database(prname, iproc, iregion_code, LOCAL_PATH)
    ! create the name of the database for the mesher and the solver
    implicit none

    integer, intent(in) :: iproc, iregion_code
    character(len=512), intent(in) :: LOCAL_PATH
    character(len=512), intent(out) :: prname

    ! name of the database file
    character(len=512) :: procname

    ! create the name for the database of the current slide and region
    write(procname,"('/proc',i6.6,'_reg',i1,'_')") iproc,iregion_code

    ! create full name with path
    prname = trim(LOCAL_PATH) // procname

  end subroutine create_name_database

  subroutine get_sys_args(kernel_parfile,kernel_file, src_mask_dir, outputfn)
    character(len=*), intent(inout) :: kernel_file, src_mask_dir,outputfn,kernel_parfile

    call getarg(1, kernel_parfile)
    call getarg(2, kernel_file)
    call getarg(3, src_mask_dir)
    call getarg(4, outputfn)

    if(trim(kernel_file) == '' .or. trim(src_mask_dir) == '' &
        .or. trim(outputfn) == '' .or. trim(kernel_parfile) == '') then
      call exit_mpi("Usage: xsrc_mask kernel_file src_mask_dir outputfn")
    endif

    if(myrank == 0) then
      write(*, *) "Kernel file  (in): ", trim(kernel_file)
      write(*, *) "Src mask dir (in): ", trim(src_mask_dir)
      write(*, *) "Output file (out): ", trim(outputfn)
    endif
  end subroutine get_sys_args

  subroutine read_event_file(eventfile, nevent, kernel_list)
    character(len=*), intent(in) :: eventfile
    integer, intent(inout) :: nevent
    character(len=*), dimension(:), allocatable, intent(inout) :: kernel_list

    ! local variables
    character(len=500) :: kernel_file
    integer :: ios, i
    integer, parameter :: fh = 1001

    open(unit=fh, file=trim(eventfile), status='old', iostat=ios)
    if ( ios /= 0 ) then
      print*, 'ERROR OPENING', trim(eventfile)
      stop
    end if

    if(myrank == 0) write(*, *) "Reading event mask: ", trim(eventfile)
    read(fh, *) nevent
    if(myrank == 0) write(*, *) "Number of events: ", nevent

    allocate(kernel_list(nevent))

    do i=1, nevent
      read(fh, '(A)') kernel_file
      if (myrank == 0) write(*, '(A, I5, A, I5, A, A)') &
        "[", i, "/", nevent, "]", trim(kernel_file)
      kernel_list(i) = trim(kernel_file)
    enddo
    close(fh)
  end subroutine read_event_file

end module src_mask_subs



program apply_src_mask
! Adios implementation: Matthieu & Ebru
! Princeton, August 2013
  use mpi
  use adios_read_mod
  use global_var, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank
  use global_var, only : init_mpi, exit_mpi, max_all_all_cr, min_all_all_cr,init_kernel_par
  use global_var, only : KERNEL_NAMES_GLOB,MODEL_NAMES_GLOB,MODEL_PERTURB_NAMES_GLOB
  use global_var, only : NPAR_GLOB,NKERNEL_GLOB,KER_HESS_NAMES_GLOB
  use AdiosIO
  use src_mask_subs
  implicit none

  ! for transverse anisotropy
  ! for isotropy

  character(len=512):: kernel_file, src_mask_dirs_file, outputfn

  ! integer, parameter :: NKERNELS = 7    !bulk_betah, bulk_betav, bulk_c, eta, hess
  ! integer, parameter :: hess_idx = 5
  ! character(len=512), parameter :: kernel_names(NKERNELS) = &
  !   (/character(len=512) :: "bulk_c_kl_crust_mantle", &
  !                           "bulk_betav_kl_crust_mantle", &
  !                           "bulk_betah_kl_crust_mantle", &
  !                           "eta_kl_crust_mantle",        &
  !                           "hess_kl_crust_mantle", &
  !                           "hess_kappa_kl_crust_mantle", &
  !                           "hess_mu_kl_crust_mantle"/)



  ! real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: kernels
  ! real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC) :: src_mask

  character(len=500), dimension(:), allocatable :: kernel_list
  integer:: nevent,ievent
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable:: kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable:: src_mask

  character(len=512) :: mask_file,src_mask_path,kernel_parfile
  real(kind=CUSTOM_REAL) :: maxv, minv
  integer :: i, ier, IIN=321

  call init_mpi()

  call get_sys_args(kernel_parfile,kernel_file, src_mask_dirs_file, outputfn)
  call init_kernel_par(kernel_parfile)

  allocate(kernels(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB*2), &
       src_mask(NGLLX, NGLLY, NGLLZ, NSPEC),stat=ier)

  if (ier .ne. 0) write(*,*) "Error in memory allocation"
  !if (kernel_names(hess_idx) /= "hess_kl_crust_mantle") call exit_mpi("Incorrect hess_idx")

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  ! read in kernel
  if (myrank == 0) write(*, *) 'Reading in kernel'
  call read_bp_file_real(kernel_file, KER_HESS_NAMES_GLOB, kernels)

  ! read in source mask
  call read_event_file(src_mask_dirs_file, nevent, kernel_list)

  do ievent=1, nevent
    src_mask_path = kernel_list(ievent)

    if (myrank==0) write(*,*) 'Reading in mask for [', ievent, &
         "/", nevent, "]: ", trim(src_mask_path)

    call create_name_database(mask_file, myrank, 1, src_mask_path)
    mask_file = trim(mask_file)//"mask_source.bin"
  if (myrank == 0) write(*, *) '[Rank0] Reading source mask file: ', trim(mask_file)
  open(IIN,file=trim(mask_file),status='old',form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(mask_file)
    call exit_mpi('file not found')
  endif
  read(IIN) src_mask
  close(IIN)

  call min_all_all_cr(minval(src_mask), minv)
  call max_all_all_cr(maxval(src_mask), maxv)
  if(myrank == 0) write(*, *) 'Min and Max value of source mask:', minv, maxv

  ! apply the source mask to kernels, excluding hessian ones
  do i=1, NKERNEL_GLOB*2
        kernels(:, :, :, :, i) = kernels(:, :, :, :, i) * src_mask
  enddo

enddo
  ! write the kernels out
  call write_bp_file(kernels, KER_HESS_NAMES_GLOB, "KERNEL_GROUPS", outputfn)

  if (myrank==0) write(*, *) 'Done applying source mask'

  call adios_finalize (myrank, ier)

  deallocate(kernels,src_mask)
  call MPI_FINALIZE(ier)

end program apply_src_mask
