subroutine get_sys_args(kernel_parfile,grad_file, precond_file, direction_file,solver_file)

  use global_var, only : myrank, exit_mpi

  character(len=512), intent(inout):: grad_file, precond_file, direction_file,solver_file,kernel_parfile

  call getarg(1, kernel_parfile)
  call getarg(2, grad_file)
  call getarg(3, solver_file)
  call getarg(4, direction_file)
  call getarg(5, precond_file)

  if(trim(grad_file) == '' .or. trim(direction_file) == '' &
       .or. trim(kernel_parfile) == '' .or. trim(solver_file) == '') then
        call exit_mpi('Usage: xsd_direction gradient_file precond_file direction_file')
  endif

  if(myrank == 0) then
    write(*, *) "Gradient file   (input): ", trim(grad_file)
    write(*, *) "Solver file  (input): ", trim(solver_file)
    write(*, *) "Direction file (output): ", trim(direction_file)
    write(*, *) "Precond file    (input): ", trim(precond_file)
  endif

end subroutine get_sys_args

program main
  use mpi
  use adios_read_mod
  use global_var
  use AdiosIO

  ! integer, parameter :: NKERNELS = 4
  ! character(len=512), parameter :: kernel_names(NKERNELS) = &
  !   (/character(len=512) :: "bulk_c_kl_crust_mantle",&
  !                           "bulk_betav_kl_crust_mantle",&
  !                           "bulk_betah_kl_crust_mantle",&
  !                           "eta_kl_crust_mantle"/)

  ! character(len=512), parameter :: precond_names(NKERNELS) = &
  !   (/character(len=512) :: "precond_bulk_c",&
  !                           "precond_bulk_betav",&
  !                           "precond_bulk_betah",&
  !                           "precond_eta"/)

 ! real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: gradient
 ! real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: precond


  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable:: gradient,precond
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable:: direction
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: jacobian
  !real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNELS):: direction

  character(len=512) :: grad_file, precond_file, direction_file,solver_file,kernel_parfile
  real(kind=CUSTOM_REAL):: gtp,gtg,gtp_old,gtg_old,gtp_all_tmp,gtg_all_tmp,max_direction,max_direction_all
  real(kind=CUSTOM_REAL),dimension(:),allocatable::gtp_all,gtg_all
  real(kind=CUSTOM_REAL)::min_direction,min_direction_all

  integer:: ier,iker

  call init_mpi()

  call get_sys_args(kernel_parfile,grad_file, precond_file, direction_file,solver_file)
  
  call init_kernel_par(kernel_parfile)

  allocate(gradient(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB),direction(NGLLX, NGLLY, NGLLZ, NSPEC,NKERNEL_GLOB))
  allocate(jacobian(NGLLX, NGLLY, NGLLZ, NSPEC))
  allocate(gtp_all(NKERNEL_GLOB),gtg_all(NKERNEL_GLOB))
  allocate(precond(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB))
  

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)

  if (myrank == 0) write(*, *) "|<----- Reading Gradient ----->| ", trim(grad_file)
  call read_bp_file_real(grad_file,KERNEL_NAMES_GLOB, gradient)

  if (myrank == 0) write(*, *) "|<----- Calculate Jacobian ----->|"
  call calculate_jacobian_matrix(solver_file, jacobian)

  if (.not. precond_file == '') then
     if (myrank == 0) write(*, *) "|<----- Reading Preconditioner ----->| ", trim(precond_file)
     call read_bp_file_real(precond_file, HESS_NAMES_GLOB, precond)
     ! steep descent method with preconditioner applied
     direction = - precond * gradient
  else
          if (myrank == 0) write(*, *) "|<----- Preconditioner = I ----->| "
     direction = -1.0 * gradient
  endif


  do iker=1,NKERNEL_GLOB
     gtg_old = sum(gradient(:,:,:,:,iker) * gradient(:,:,:,:,iker))
     gtp_old = sum(direction(:,:,:,:,iker) * gradient(:,:,:,:,iker))
     call mpi_barrier(MPI_COMM_WORLD, ier)
     call sum_all_all_cr(gtg_old, gtg_all_tmp)
     call sum_all_all_cr(gtp_old, gtp_all_tmp)
     gtg_all(iker) = gtg_all_tmp
     gtp_all(iker) = gtp_all_tmp
  enddo

  gtp_old = sum(gtp_all)
  gtg_old = sum(gtg_all)


  do iker = 1,NKERNEL_GLOB
     call max_all_all_cr(maxval(direction(:, :, :, :,iker)), max_direction)
     call min_all_all_cr(minval(direction(:, :, :, :,iker)), min_direction)
     if (myrank == 0) then                                                 
        write(*,*) " Gradient Stats (After Preconditioner)"                
        write(*, '(A, ES12.2, 5X, ES12.2)') &                              
             trim(KERNEL_NAMES_GLOB(iker))//"(max,min)", &                 
             max_direction, min_direction                                  
     endif
  enddo


  max_direction = maxval(abs(direction))
  call max_all_all_cr(max_direction,max_direction_all)

  min_direction = minval(abs(direction))
  call min_all_all_cr(min_direction,min_direction_all)


  
  call Parallel_ComputeInnerProduct(gradient, gradient, &
       NKERNEL_GLOB, jacobian, gtg)

  call Parallel_ComputeInnerProduct(gradient, direction, &
       NKERNEL_GLOB, jacobian, gtp)

  
  open(99, file = "gtg")  
  write(99,*) gtg_old
  close(99)

  open(99, file = "gtp")  
  write(99,*) gtp_old
  close(99)

  !#######
  if (myrank == 0) write(*, *) "|<----- Writing Direction ----->| ", trim(direction_file)
  call write_bp_file(direction, KERNEL_NAMES_GLOB, "KERNEL_GROUPS", direction_file)

  call adios_finalize(myrank, ier)

  deallocate(gradient,direction)
  deallocate(jacobian,gtp_all,gtg_all)

  if (myrank==0) then
  write(*,*) "gtp_old : ",gtp_old," gtg_old :",gtg_old
  write(*,*) "gtp : ",gtp," gtg :",gtg
  write(*,*) "min_val_direction : " , min_direction_all," max_val_direction: ", max_direction_all
endif

  call MPI_FINALIZE(ier)

end program
