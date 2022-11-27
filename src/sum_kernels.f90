module sum_kernels_subs
  use mpi
  use global_var, only : CUSTOM_REAL, myrank, exit_mpi,init_kernel_par,nprocs,CUSTOM_MPI_TYPE
  use global_var, only : NGLLX, NGLLY, NGLLZ, NSPEC,NPAR_GLOB,ATTENUATION_FLAG,NKERNEL_GLOB
  implicit none
   ! model updates
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: models

contains


  subroutine init_sum_kernel_model()
    integer::ier
    ! Allocates memory for variables in model_par
    allocate(models(NGLLX,NGLLY,NGLLZ,NSPEC,NPAR_GLOB),stat=ier)

    models=0.0d0
  
    if (ier /= 0) stop " Error Allocating parameter for Sum Kernel : Init Model Variables "

  end subroutine init_sum_kernel_model

  subroutine clean_sum_kernel_model()
    deallocate(models)
    
  end subroutine clean_sum_kernel_model

  subroutine read_event_file(eventfile, nevent, kernel_list, weight_list)
    character(len=*), intent(in) :: eventfile
    integer, intent(inout) :: nevent
    character(len=*), dimension(:), allocatable, intent(inout) :: kernel_list
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: weight_list

    ! local variables
    character(len=500) :: kernel_file
    real(kind=CUSTOM_REAL) :: weight
    integer :: ios, i
    integer, parameter :: fh = 1001

    open(unit=fh, file=trim(eventfile), status='old', iostat=ios)
    if ( ios /= 0 ) then
      print*, 'ERROR OPENING', trim(eventfile)
      stop
    end if

    if(myrank == 0) write(*, *) "Reading events and weights from: ", trim(eventfile)
    read(fh, *) nevent
    if(myrank == 0) write(*, *) "Number of events: ", nevent

    allocate(kernel_list(nevent))
    allocate(weight_list(nevent))

    do i=1, nevent
      read(fh, *) weight
      read(fh, '(A)') kernel_file
      if (myrank == 0) write(*, '(A, I5, A, I5, A, A, A, ES16.8)') &
        "[", i, "/", nevent, "]", trim(kernel_file), ": ", weight
      kernel_list(i) = trim(kernel_file)
      weight_list(i) = weight
    enddo
    close(fh)
  end subroutine read_event_file

  subroutine get_sys_args(kernel_parfile,eventfile, outputfn,inputmodel,eventfile_qmu)
    character(len=*), intent(inout) :: eventfile, outputfn,inputmodel,kernel_parfile
    character(len=*), intent(inout) :: eventfile_qmu
    call getarg(1, kernel_parfile)
    call getarg(2, eventfile)
    call getarg(3, outputfn)
    call getarg(4,inputmodel)
    call getarg(5,eventfile_qmu)
    
    if(trim(eventfile) == '' .or. trim(outputfn) == '' .or. trim(kernel_parfile) == '') then
      call exit_mpi("Usage: xsum_kernels eventfile outputfn")
   endif

    if(myrank == 0) then
      write(*, *) "Event file(in): ", trim(eventfile)
      write(*, *) "Output file(out, kernel sums): ", trim(outputfn)
      if (len(trim(eventfile_qmu)) .gt. 0) write(*, *) "Event file attenuation(in): ", trim(eventfile_qmu) 
      write(*, *) "Input model(in): ", trim(inputmodel)
    endif
  end subroutine get_sys_args


    subroutine clip_sem(kernels,percentile,nvar)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: kernels
    real(kind=CUSTOM_REAL),intent(in) :: percentile
    integer,intent(in) :: nvar
    real(kind=CUSTOM_REAL), dimension(:,:),allocatable:: ksort,ksort_full
    real(kind=CUSTOM_REAL),dimension(:),allocatable::kpercent
    real(kind=CUSTOM_REAL) :: val
    integer :: ier,i,j,k,ispec,iker,npercent


    npercent = int(floor(percentile * NSPEC * nprocs))
    
    allocate(ksort(NSPEC,nvar),stat=ier)
    allocate(ksort_full(NSPEC*nprocs,nvar),kpercent(nvar),stat=ier)

    ksort = 0._CUSTOM_REAL
    ksort_full = 0._CUSTOM_REAL
    kpercent = 0._CUSTOM_REAL
    ! Prepare array for sorting
    do iker = 1,nvar
       do ispec = 1, NSPEC 
          ksort(ispec,iker) = maxval(abs(kernels(:,:,:,ispec,iker)))
       enddo
       call quicksort(ksort(:,iker),1,NSPEC)
       call MPI_ALLGATHER(ksort(:,iker),NSPEC,CUSTOM_MPI_TYPE,ksort_full(:,iker),NSPEC, &
            CUSTOM_MPI_TYPE,MPI_COMM_WORLD,ier)
       call MPI_Barrier(MPI_COMM_WORLD, ier)
       call quicksort(ksort_full(:,iker),1,NSPEC*nprocs)
       kpercent(iker) = ksort_full(npercent,iker)
    enddo

    if (myrank==0) then
       write(*,*) "Percentil Kernels"
       do i=1,nvar
          write(*,*) "<Percentile,Value Kernel>",percentile*100,kpercent(i)
       enddo
    endif

    do iker = 1,nvar
       do ispec = 1, NSPEC
          do k = 1, NGLLZ
             do j = 1, NGLLY
                do i = 1, NGLLX
                   val = abs(kernels(i,j,k,ispec,iker))
                   if (val > kpercent(iker)) then
                      kernels(i,j,k,ispec,iker) = 0._CUSTOM_REAL
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
 
    call MPI_Barrier(MPI_COMM_WORLD, ier)
    
  end subroutine clip_sem


  recursive subroutine quicksort(a, first, last)
    implicit none
    real(kind=CUSTOM_REAL),dimension(:),intent(inout) :: a
    integer,intent(in) :: first, last
    real :: x, t
    integer :: i, j

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do
    if (first < i-1) call quicksort(a, first, i-1)
    if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort
  

end module sum_kernels_subs


program sum_kernels
  use mpi
  use adios_read_mod
  use global_var, only : CUSTOM_REAL, NGLLX, NGLLY, NGLLZ, NSPEC, myrank,init_kernel_par,NPAR_GLOB,NKERNEL_GLOB
  use global_var, only : init_mpi, exit_mpi,KERNEL_NAMES_GLOB,MODEL_NAMES_GLOB,KERNEL_NAMES_GLOB_NQ,NHESS0
  use global_var, only : ATTENUATION_FLAG,QMU_IDX,KQMU_IDX,KER_HESS_NAMES_GLOB,max_all_all_cr,HESS_NAMES_GLOB
  use AdiosIO
  use sum_kernels_subs

  implicit none

  integer:: nevent, ievent, ier,i,idx
  character(len=500) :: eventfile, outputfn, kernel_file,input_model_file,kernel_parfile,kernel_file_qmu
  character(len=500) :: eventfile_qmu
  character(len=500), dimension(:), allocatable :: kernel_list,kernel_list_qmu
  real(kind=CUSTOM_REAL):: weight
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: weight_list,weight_list_qmu

  !integer, parameter :: NKERNELS = 5    !bulk_betah, bulk_betav, bulk_c, eta, rho, hess
  !integer, parameter :: hess_idx = 5
  !character(len=500), parameter :: kernel_names(NKERNELS) = &
  !  (/character(len=150) :: "bulk_betah_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
  !                          "bulk_c_kl_crust_mantle",     "eta_kl_crust_mantle",        &
  !                          "hess_kl_crust_mantle"/)

  ! integer, parameter :: NKERNELS = 7    !bulk_betah, bulk_betav, bulk_c, eta, rho, hess
  ! integer, parameter :: hess_idx = 5
  ! character(len=500), parameter :: kernel_names(NKERNELS) = &
  !   (/character(len=500) :: "bulk_betah_kl_crust_mantle", "bulk_betav_kl_crust_mantle", &
  !                           "bulk_c_kl_crust_mantle", "eta_kl_crust_mantle",        &
  !                           "hess_kl_crust_mantle", "hess_kappa_kl_crust_mantle", &
  !                           "hess_mu_kl_crust_mantle"/)

  ! real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: total_kernel
  ! real(kind=CUSTOM_REAL),dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: kernels


  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: total_kernel
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:),allocatable:: kernels,hessian

  call init_mpi()

  call get_sys_args(kernel_parfile,eventfile, outputfn,input_model_file,eventfile_qmu)
  
  call init_kernel_par(kernel_parfile)

  allocate(total_kernel(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB*2),kernels(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB))
  allocate(hessian(NGLLX, NGLLY, NGLLZ, NSPEC, NHESS0))


!  if (kernel_names(hess_idx) /= "hess_kl_crust_mantle") call exit_mpi("Incorrect hess_idx")

  
  call read_event_file(eventfile, nevent, kernel_list, weight_list)
  
      

  if ((ATTENUATION_FLAG) .and. (len(trim(eventfile_qmu)) .gt. 0)) then
     if (myrank==0) write(*,*) " ---> Attenuation event file <---"
     call read_event_file(eventfile_qmu, nevent, kernel_list_qmu, weight_list_qmu)
  end if
  

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
       "verbose=1", ier)


  total_kernel = 0.
  do ievent=1, nevent
     
     kernel_file = kernel_list(ievent)
     weight = weight_list(ievent)


     if ((ATTENUATION_FLAG) .and. (NPAR_GLOB .gt. 1)) then
     !   if (myrank ==0) then
     !       do i=1,NKERNEL_GLOB -1
     !         print*, "reading: ",trim(KERNEL_NAMES_GLOB_NQ(i)) 
     !       enddo
     !    endif

        if (myrank==0) write(*,*) 'Reading in kernel for [', ievent, &
             "/", nevent, "]: ", trim(kernel_file), " | weight: ", weight
             
        call read_bp_file_real(kernel_file, KERNEL_NAMES_GLOB_NQ, kernels)
     else
        call read_bp_file_real(kernel_file, KERNEL_NAMES_GLOB, kernels)
     endif
     

    if ((ATTENUATION_FLAG) .and. (len(trim(eventfile_qmu)) .gt. 0)) then
       kernel_file_qmu = kernel_list_qmu(ievent)
       weight = weight_list_qmu(ievent)
       if (myrank==0) write(*,*) 'Reading in kernel (Attenuation) for [', ievent, &
              "/", nevent, "]: ", trim(kernel_file_qmu), " | weight: ", weight
       call read_bp_file_real(kernel_file_qmu,KERNEL_NAMES_GLOB(KQMU_IDX),kernels(:,:,:,:,KQMU_IDX))
       
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ier)
    ! only keep the abs value of the hess kernel
    !kernels(:, :, :, :, hess_idx) = abs(kernels(:, :, :, :, hess_idx))
    !kernels(:, :, :, :, hess_idx:(hess_idx+2)) = abs(kernels(:, :, :, :, hess_idx:(hess_idx+2)))

    if (myrank == 0) write(*, *) 'Cliping Kernels 99.8',ievent
    call clip_sem(kernels,0.9950,NKERNEL_GLOB)
    
    do idx=1,NKERNEL_GLOB
       total_kernel(:,:,:,:,idx) = total_kernel(:,:,:,:,idx) + kernels(:,:,:,:,idx) * weight
    enddo

    if (NHESS0 > 0 ) then
       do idx=1,NHESS0
          if (HESS_NAMES_GLOB(idx)(1:4) .ne. &
               "_qmu") then
             call read_bp_file_real(kernel_file, &
                  HESS_NAMES_GLOB(idx)(5:len_trim(HESS_NAMES_GLOB(idx))), &
                  hessian(:,:,:,:,idx))
          endif
       enddo
       
       if (myrank == 0) write(*, *) 'Cliping Hessian 99.8',ievent
       call clip_sem(hessian,0.9980,NKERNEL_GLOB)
     
       do idx=NKERNEL_GLOB+1,NKERNEL_GLOB + NHESS0
          total_kernel(:,:,:,:,idx) = &
               total_kernel(:,:,:,:,idx) + &
               abs(hessian(:,:,:,:,idx - NKERNEL_GLOB)) * weight 
       enddo
             

       if ((ATTENUATION_FLAG) .and. (len(trim(eventfile_qmu)) .gt. 0)) then
          if (myrank==0) write(*,*) 'Reading in Hessian (Attenuation)'
          if (ievent==1) total_kernel(:,:,:,:,KQMU_IDX + NKERNEL_GLOB) = 0.0d0
          call read_bp_file_real(kernel_file_qmu, &
               HESS_NAMES_GLOB(KQMU_IDX)(5:len_trim(HESS_NAMES_GLOB(KQMU_IDX))), &
               hessian(:,:,:,:,KQMU_IDX))
          total_kernel(:,:,:,:,KQMU_IDX + NKERNEL_GLOB)  = &
               total_kernel(:,:,:,:,KQMU_IDX + NKERNEL_GLOB) + &
               abs(hessian(:,:,:,:,KQMU_IDX)) * weight
       endif
    endif
    
 enddo


 call write_bp_file(total_kernel, KER_HESS_NAMES_GLOB, "KERNEL_GROUPS", outputfn)
 

  if (ATTENUATION_FLAG) then
     if (myrank==0) print*, "Reading model to cast dq to logdq"
     call init_sum_kernel_model()
     ! reads in current transverse isotropic model files: vpv.. & vsv.. & eta & rho
     if(myrank == 0) print*, "|<---- Reading Model File ---->|", trim(input_model_file)
     call read_bp_file_real(input_model_file, MODEL_NAMES_GLOB, models)

     total_kernel = 0.0d0
     
     if (myrank == 0) write(*, *) "|<----- Reading Sum Kernel ----->| ", trim(outputfn)
     call read_bp_file_real(outputfn,KER_HESS_NAMES_GLOB, total_kernel)

     total_kernel(:,:,:,:,KQMU_IDX) = &
          total_kernel(:,:,:,:,KQMU_IDX) * 1.0 / models(:,:,:,:,QMU_IDX)
     
     ! total_kernel(:,:,:,:,KQMU_IDX + NKERNEL_GLOB) = &
     !      total_kernel(:,:,:,:,KQMU_IDX + NKERNEL_GLOB) &
     !      * (1.0 / models(:,:,:,:,QMU_IDX))**2.0


     call write_bp_file(total_kernel, KER_HESS_NAMES_GLOB, "KERNEL_GROUPS", outputfn)
     
     

     call clean_sum_kernel_model()
     
  endif

  

  if (myrank==0) print*, 'Done summing all the kernels'
  if (myrank==0) print*, "output file: ", trim(outputfn)

  call MPI_Barrier(MPI_COMM_WORLD, ier)
  
  call adios_finalize(myrank, ier)
  
  deallocate(total_kernel,kernels)

  call MPI_Barrier(MPI_COMM_WORLD, ier)
  call MPI_FINALIZE(ier)

end program sum_kernels
