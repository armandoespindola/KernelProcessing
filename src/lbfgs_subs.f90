module lbfgs_subs

  use mpi
  use global_var, only : CUSTOM_REAL, myrank, exit_mpi, NGLLX, NGLLY, NGLLZ, &
                     NSPEC,NPAR_GLOB,KERNEL_NAMES_GLOB,NSPEC,NKERNEL_GLOB,MODEL_PERTURB_NAMES_GLOB
  use global_var, only : Parallel_ComputeInnerProduct, Parallel_ComputeL2normSquare,sum_all_all_cr
  use global_var, only : min_all_all_cr,max_all_all_cr
  use AdiosIO, only : read_bp_file_real
  implicit none

  contains

  subroutine get_sys_args(kernel_parfile,input_path_file, solver_file,precond_file, outputfn)
    character(len=*), intent(inout) :: input_path_file, solver_file, outputfn,kernel_parfile
    character(len=*), intent(inout) :: precond_file

    if(myrank == 0) print*, '|<============= Get System Args =============>|'

    call getarg(1, kernel_parfile)
    call getarg(2, input_path_file)
    call getarg(3, solver_file)
    call getarg(4, outputfn)
    call getarg(5, precond_file)

    if(trim(input_path_file) == '' .or. trim(solver_file) == '' &
        .or. trim(outputfn) == '' .or. trim(kernel_parfile) == '') then
      call exit_mpi("Usage: ./xlbfgs input_path_file solver_file outputfn H0_file[optional]")
    endif

    if(myrank == 0) then
      print*, "Input path file: ", trim(input_path_file)
      print*, "Solver file:", trim(solver_file)
      print*, "Precond file: ", trim(precond_file)
      print*, "Output file: ", trim(outputfn)
      
    endif
  end subroutine get_sys_args



  subroutine get_sys_args_bm(kernel_parfile,input_path_file, solver_file,perturb_file,precond_file, outputfn)
    character(len=*), intent(inout) :: input_path_file, solver_file, outputfn,kernel_parfile
    character(len=*), intent(inout) :: precond_file,perturb_file

    if(myrank == 0) print*, '|<============= Get System Args =============>|'

    call getarg(1, kernel_parfile)
    call getarg(2, input_path_file)
    call getarg(3, solver_file)
    call getarg(4, outputfn)
    call getarg(5, perturb_file)
    call getarg(6, precond_file)

    if(trim(input_path_file) == '' .or. trim(solver_file) == '' &
         .or. trim(outputfn) == '' .or. trim(kernel_parfile) == '' &
         .or. trim(perturb_file) == '') then
      call exit_mpi("Usage: ./xlbfgs_bm input_path_file solver_file outputfn H0_file[optional]")
    endif

    if(myrank == 0) then
      print*, "Input path file: ", trim(input_path_file)
      print*, "Solver file:", trim(solver_file)
      print*, "Precond file: ", trim(precond_file)
      print*, "Output file: ", trim(outputfn)
      print*, "Perturb file: ", trim(perturb_file)
      
    endif
  end subroutine get_sys_args_bm

  subroutine parse_input_path(fn, niter, gradient_files, model_change_files)
    character(len=*), intent(in) :: fn
    integer, intent(inout) :: niter
    character(len=512), dimension(:), allocatable, intent(inout) :: gradient_files, &
                                                                    model_change_files

    integer :: fh = 1001
    integer :: i, ier

    if(myrank == 0) print*, '|<============= Parsing Input Path File =============>|'

    open(fh, file=trim(fn), status='old')

    read(fh, *) niter
    if(myrank == 0) print*, "niter: ", niter
    ! The gradient files is for each iteration, gradient of M(n), ..., M(n-k)
    ! So the size is (k + 1)
    allocate(gradient_files(niter+1))
    ! The model change files is the model change between two continues
    ! So the size is (k)
    allocate(model_change_files(niter))

    ! read the last k iteration
    do i=1, niter
      read(fh, '(A)') gradient_files(i)
      read(fh, '(A)') model_change_files(i)
      if(myrank == 0) then
        write(*, '(A, I3, A, A)'), "[iter=", i, "] grad file:   ", trim(gradient_files(i))
        write(*, '(A, I3, A, A)'), "[iter=", i, "] dkernel file:", trim(model_change_files(i))
      endif
    enddo

    ! read the gradient from current iteration
    read(fh, '(A)') gradient_files(niter+1)
    if(myrank == 0) then
      print*, "[Current iter] grad file: ", trim(gradient_files(niter+1))
    endif

    close(fh)

    call MPI_Barrier(MPI_COMM_WORLD, ier)
  end subroutine parse_input_path

  subroutine read_all_bp_files(niter, nkernels, gradient_files, model_change_files, &
                               kernel_names, gradient, yks, sks)

    integer, intent(in) :: niter
    character(len=*), dimension(:), intent(in) :: gradient_files, model_change_files
    integer, intent(in) :: nkernels
    character(len=*), dimension(:), intent(in) :: kernel_names

    ! (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :) :: gradient

    ! (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels, niter)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :) :: yks, sks

    ! local variable
    ! (NGLLX, NGLLY, NGLLZ, NSPEC, nkernels, niter)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: stored_gradients

    integer :: i
    
    !character(len=500), dimension(NKERNEL_GLOB) :: names_dmodel

    !names_dmodel=MODEL_PERTURB_NAMES_GLOB(QMU_IDX)


    if(myrank == 0) print*, '|<============= Reading Previous BP files =============>|'
    allocate(stored_gradients(NGLLX, NGLLY, NGLLZ, NSPEC, nkernels, niter+1))

    ! Read gradient
    do i=1, niter+1
      if(myrank == 0) write(*, '(A, I2, A, A)') &
        "Reading [iter=", i, "] Gradient: ", trim(gradient_files(i))
      call read_bp_file_real(gradient_files(i),KERNEL_NAMES_GLOB, &
                             stored_gradients(:, :, :, :, :, i))
    enddo

    ! store the last iteration gradient
    gradient = stored_gradients(:, :, :, :, :, niter+1)

    ! Calculate the gradient change
    do i=1, niter
      yks(:, :, :, :, :, i) = &
        stored_gradients(:, :, :, :, :, i+1) - stored_gradients(:, :, :, :, :, i)
    enddo

    ! Read Model Change
    do i=1, niter
      if(myrank == 0) write(*, '(A, I2, A, A)') &
        "Reading [iter=", i, "] dkernel: ", trim(model_change_files(i))
      call read_bp_file_real(model_change_files(i), KERNEL_NAMES_GLOB, &
                             sks(:, :, :, :, :, i))
    enddo

    deallocate(stored_gradients)

  end subroutine read_all_bp_files

  subroutine calculate_LBFGS_direction(niter, nkernels, jacobian, gradient, &
                                       precond, yks, sks, direction)
    integer, intent(in) :: niter
    integer, intent(in) :: nkernels
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: gradient
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: precond
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), intent(in) :: yks, sks
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: direction

    integer :: i
    real(kind=CUSTOM_REAL) :: tmp, beta
    real(kind=CUSTOM_REAL) :: norm_y,rhok
    real(kind=CUSTOM_REAL), dimension(:), allocatable :: pk_store, ak_store

    if(myrank == 0) print*, '|<============= Calculating LBFGS Direction =============>|'
    allocate(pk_store(niter))
    allocate(ak_store(niter))

    direction = gradient

    ! first round
    do i=niter, 1, -1
      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
                                        yks(:, :, :, :, :, i), &
                                        nkernels, jacobian, tmp)
      pk_store(i) = 1.0 / tmp

      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
                                        direction, &
                                        nkernels, jacobian, tmp)
      ak_store(i) = pk_store(i) * tmp

      if(myrank == 0) write(*, '(A, I2, A, i2, A, ES18.8, A, ES18.8)') &
        "First Loop[iter=", i, "/", niter, "] pk=", pk_store(i), ", ak=", ak_store(i)

      direction = direction - ak_store(i) * yks(:, :, :, :, :, i)
    enddo

    ! Precondition 1
    call Parallel_ComputeL2normSquare(yks(:, :, :, :, :, niter), nkernels, &
         jacobian, norm_y)

    
    rhok = 1.0 / (pk_store(niter) * norm_y)
    direction = rhok * direction * precond

    ! Precondition 2
    !direction = precond * direction

    ! second round
    do i=1, niter, 1
      call Parallel_ComputeInnerProduct(yks(:, :, :, :, :, i),  direction, &
                                        nkernels, jacobian, tmp)
      beta = pk_store(i) * tmp

      if(myrank == 0) write(*, '(A, I3, A, I3, A, ES18.8)') &
        "Second Loop[iter=", i, "/", niter, "] beta=", beta
      direction = direction + (ak_store(i) - beta) * sks(:, :, :, :, :, i)
    enddo

    ! reverse to get descent direction
    direction = -1.0 * direction

  end subroutine calculate_LBFGS_direction


    subroutine calculate_LBFGS_Bm(niter, nkernels, jacobian, test_m, &
                                       precond, yks, sks,Bm)
    integer, intent(in) :: niter
    integer, intent(in) :: nkernels
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: test_m
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: precond
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), intent(in) :: yks, sks
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: Bm

    integer :: i,j
    real(kind=CUSTOM_REAL) :: tmp1,tmp2, beta,gamma,sigma,bts,ats,tmp,norm_y
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), allocatable :: ak,bk

    if(myrank == 0) print*, '|<============= Calculating LBFGS Bm =============>|'
    allocate(ak(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB, niter))
    allocate(bk(NGLLX, NGLLY, NGLLZ, NSPEC, NKERNEL_GLOB, niter))

    ak = 0.0d0
    bk = 0.0d0
    Bm = 0.0d0



    call Parallel_ComputeL2normSquare(sks(:, :, :, :, :, niter), nkernels, &
         jacobian, norm_y)

    call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, niter), &
                                        yks(:, :, :, :, :, niter), &
                                        nkernels, jacobian, tmp)


    gamma = tmp / norm_y
    sigma = gamma
    !sigma = 1.0d0 / gamma

    ! first round
    do i=1, niter
       
      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
                                        yks(:, :, :, :, :, i), &
                                        nkernels, jacobian, tmp)

      bk(:, :, :, :, :, i) = yks(:, :, :, :, :, i) / sqrt(tmp)



      ! B0 approximation
      ak(:, :, :, :, :, i) = sks(:, :, :, :, :, i) * sigma


      do j=1,i - 1,1
         
         call Parallel_ComputeInnerProduct(bk(:, :, :, :, :, j), &
              sks(:, :, :, :, :, i), &
              nkernels, jacobian, bts)


         call Parallel_ComputeInnerProduct(ak(:, :, :, :, :, j), &
              sks(:, :, :, :, :, i), &
              nkernels, jacobian, ats)


         ak(:, :, :, :, :, i) = ak(:, :, :, :, :, i) + (bts * bk(:, :, :, :, :, j) - &
              ats * ak(:, :, :, :, :, j))

         if(myrank == 0) write(*, '(A, I2, A, i2, A)') &
              "Interior Loop [iter=", j, "/", i-1, "]"
         
      enddo
      

      call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
                                        ak(:, :, :, :, :, i), &
                                        nkernels, jacobian, tmp)


      ak(:, :, :, :, :, i) = ak(:, :, :, :, :, i) / sqrt(tmp)
     

      if(myrank == 0) write(*, '(A, I2, A, i2, A)') &
        "Loop Approximation[iter=", i, "/", niter, "]"

   enddo

   do i=1,niter
      
         call Parallel_ComputeInnerProduct(test_m(:, :, :, :, :), &
                                        ak(:, :, :, :, :, i), &
                                        nkernels, jacobian, tmp1)


         call Parallel_ComputeInnerProduct(test_m(:, :, :, :, :), &
                                        bk(:, :, :, :, :, i), &
                                        nkernels, jacobian, tmp2)


      Bm(:, :, :, :, :) = Bm(:, :, :, :, :) + (bk(:, :, :, :, :, i) * tmp2 - &
           ak(:, :, :, :, :, i) * tmp1)


      if(myrank == 0) write(*, '(A, I2, A, i2, A)') &
        "Applying B to m (Bm)[iter=", i, "/", niter, "]"

   enddo
   

   Bm(:, :, :, :, :) = Bm(:, :, :, :, :) + sigma * test_m(:, :, :, :, :)  



  end subroutine calculate_LBFGS_Bm



    subroutine calculate_diag_SR1(niter, nkernels, jacobian, &
                                       precond, yks, sks,Hdiag)
    integer, intent(in) :: niter
    integer, intent(in) :: nkernels
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: precond
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :, :), intent(in) :: yks, sks
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: Hdiag

    integer :: i
    real(kind=CUSTOM_REAL) :: tmp, tmp2,norm_sk,norm_yk,norm_Bs,norm_Hy,norm_y,norm_sk_yk
    real(kind=CUSTOM_REAL) :: testH,testB

    if(myrank == 0) print*, '|<============= Calculating diag(H) SR1 Direction =============>|'
    call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, niter), &
                                        yks(:, :, :, :, :, niter), &
                                        nkernels, jacobian, norm_sk_yk)
    
    call Parallel_ComputeL2normSquare(yks(:, :, :, :, :, niter), nkernels, &
         jacobian, norm_y)
    
    Hdiag = precond * norm_sk_yk / norm_y

    if (myrank == 0) write(*,'(A,ES18.8)') " [ Multiple of Beta *  I ] Beta = ", &
         norm_sk_yk / norm_y

    ! first round
    do i=1, niter

       call Parallel_ComputeInnerProduct(yks(:, :, :, :, :, i), &
                                   sks(:, :, :, :, :, i) - yks(:, :, :, :, :, i) * Hdiag, &
                                   nkernels, jacobian, tmp)

       call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
            yks(:, :, :, :, :, i) - sks(:, :, :, :, :, i) * (1.0 / Hdiag), &
            nkernels, jacobian, tmp2)


       call Parallel_ComputeL2normSquare(yks(:, :, :, :, :, i), nkernels, &
            jacobian, norm_yk)

       call Parallel_ComputeInnerProduct(sks(:, :, :, :, :, i), &
                                        yks(:, :, :, :, :, i), &
                                        nkernels, jacobian, norm_sk_yk)


       call Parallel_ComputeL2normSquare(sks(:, :, :, :, :, i), nkernels, &
            jacobian, norm_sk)


       call Parallel_ComputeL2normSquare(yks(:, :, :, :, :, i) - sks(:, :, :, :, :, i) & 
            * (1.0 / Hdiag), nkernels, &
            jacobian, norm_Bs)
       
       call Parallel_ComputeL2normSquare(sks(:, :, :, :, :, i) - yks(:, :, :, :, :, i) &
            * Hdiag, nkernels, &
            jacobian, norm_Hy)


       testH = sqrt(norm_Hy * norm_yk) * 1.0E-8
       testB = sqrt(norm_Bs * norm_sk) * 1.0E-8


       if(myrank == 0) write(*, '(A, ES18.8, ES18.8, ES18.8)') &
            "Values [yTs, y^2, gamma] : ",norm_sk_yk,norm_yk,norm_sk_yk/norm_yk
       
       if(myrank == 0) write(*, '(A, I2, A, i2, A, ES18.8, ES18.8, ES18.8, A, ES18.8, ES18.8)') &
            "First Loop[iter=", i, "/", niter, "] [H] r(k) = ", abs(tmp),tmp,testH, &
            " [---] [B] r(k) =",abs(tmp2),testB 

       if (abs(tmp) > testH) then
          Hdiag = Hdiag + (sks(:, :, :, :, :, i) - yks(:, :, :, :, :, i) * Hdiag) * &
               (sks(:, :, :, :, :, i) - yks(:, :, :, :, :, i) * Hdiag) * (1.0 / tmp)
       end if
       

       ! Hdiag = Hdiag + (sks(:, :, :, :, :, i) * yks(:, :, :, :, :, i)) / &
       ! (yks(:, :, :, :, :, i)  * yks(:, :, :, :, :, i))
       

    enddo
    
  end subroutine calculate_diag_SR1



  subroutine check_status(grad,dir,jacobian,nkernels)
    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(in) :: grad
    integer, intent(in) :: nkernels
    real(kind=CUSTOM_REAL), dimension(:, :, :, :), intent(in) :: jacobian

    real(kind=CUSTOM_REAL), dimension(:, :, :, :, :), intent(inout) :: dir
    real(kind=CUSTOM_REAL) :: gg,dd,gd,angle,flag,flag_all
    real(kind=CUSTOM_REAL) :: max_direction,max_direction_all,min_direction,min_direction_all
    integer:: ier

    call Parallel_ComputeInnerProduct(grad, &
                                        grad, &
                                        nkernels, jacobian, gg)

    call Parallel_ComputeInnerProduct(dir, &
                                        dir, &
                                        nkernels, jacobian, dd)

    call Parallel_ComputeInnerProduct(grad, &
                                        dir, &
                                        nkernels, jacobian, gd)

    call mpi_barrier(MPI_COMM_WORLD, ier)
    
    if (myrank == 0) then 
    open(99, file = "gtg")  
    write(99,*) gg
    close(99)

    open(99, file = "gtp")  
    write(99,*) gd
    close(99)
    end if

    angle = acos(-1.0 * gd / (gg * dd)**0.5) * 180.0 / 3.14156d0

    if (angle .gt. 90.0d0) then
       flag=1.0
    end if

    call sum_all_all_cr(flag,flag_all)
    


    max_direction = maxval(abs(dir))
    call max_all_all_cr(max_direction,max_direction_all)

    min_direction = minval(abs(dir))
    call min_all_all_cr(min_direction,min_direction_all)

    if(myrank == 0) print*, '|<============= Check Status LBFGS =============>|'
    
    if (myrank == 0) write(*,*) "Angle search direction from current: ", angle
    
    if (flag_all .gt. 0) then 
       if(myrank == 0) print*, '|<============= Restarting LBFGS =============>|'
       dir = -1.0 * grad
       
       call Parallel_ComputeInnerProduct(grad, &
                                        grad, &
                                        nkernels, jacobian, gg)

       call Parallel_ComputeInnerProduct(dir, &
            dir, &
            nkernels, jacobian, dd)

       call Parallel_ComputeInnerProduct(grad, &
            dir, &
            nkernels, jacobian, gd)

       call mpi_barrier(MPI_COMM_WORLD, ier)
    
       if (myrank == 0) then 
          open(99, file = "gtg")  
          write(99,*) gg
          close(99)

          open(99, file = "gtp")  
          write(99,*) gd
          close(99)
       end if

       max_direction = maxval(abs(dir))
       call max_all_all_cr(max_direction,max_direction_all)

       min_direction = minval(abs(dir))
       call min_all_all_cr(min_direction,min_direction_all)
    end if



    if(myrank == 0) then 
       write(*,*) "gtp_old : ",gd," gtg_old :",gg
       write(*,*) "min_val_direction : " , min_direction_all," max_val_direction: ", max_direction_all
    end if

   

  end subroutine check_status

end module lbfgs_subs
