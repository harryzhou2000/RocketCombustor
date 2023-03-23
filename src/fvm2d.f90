!!==================================================================
!!  Parallel Finite Volume Solver for 2D Rocket Combustion,
!!  Qian Wang, 2020/05
!!  Email: qian.wang@epfl.ch
!!------------------------------------------------------------------

program main
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer i,k,rk,istep,istop,inner_rk,ierr
    real*8 :: totaltime_num,time_all,time0,time1,ttime_n,RoRes,RoRes_tmp

    call startup

    call ReadMesh
    
    write(*,*) 'Flowfield initialization begin at processor',myProc

    call initialization
    if(IfReadBackup) call Read_Backup
    
    call compute_reconstruction_matrix        !! second or high order reconstruction initialization

    call detect_common_boundary

    !! for time-dependent boundary condition
    if(IfRightDomain) then
      call arrange_common_boundary_cell_index
      call load_neural_network
    endif

    !! output init flow-field
    call output_vtk(0)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    write(*,*) 'Reconstruction initialization finished at processor', myProc
    
    write(*,*) 'begin computation loop at processor', MyProc

    totaltime_num=totaltime*(1.-1.e-10)
    
    istop=0;  time_all=0.

    call cpu_time(time0)
    
    do istep=1, nstep_max

      if(IfExplicit) call timestep

      if( ttime+dt> totaltime_num )then
        dt=totaltime-ttime
        if(IfExplicit) step(:)=dt
        istop=1
      endif

      ! write(*,*) 'time step computation finished at processor', MyProc

      ttime_n=ttime
      
      do rk=1, InLoop
      
        if(.not. IfSteady) then
          ttime= ttime_n + cs(rk)*dt   !! for explicit, cs(rk-1); for implicit, cs(rk) 
          call compute_p_out
          if(IfRightDomain) call compute_common_boundary_value
        endif

        !!! steady-state problem
        if(IfSteady) then  

          if(IfExplicit)then
            call erk_stage(istep,rk)  !! third order TVD Runge-Kutta method / RK4, local time stepping
          else
            call lu_sgs_stage(istep)  !! implicit LU-SGS
          endif
        
        !!! unsteady problem
        else               

          if(IfExplicit)then
            call erk_stage(istep,rk)  !! third order TVD Runge-Kutta method / RK4
          else
            select case(implicit_time_method)
            case(1)        ! implicit RK
              call irk_stage(istep,rk)
            case(2)        ! trapezoidal rule 
              call trapezoidal_stage(istep)
            case(3)        ! backward Euler 
              call backwardEuler_stage(istep)
            case(4)        ! implicit midpoint 
              call implicitMidpoint_stage(istep)
            case(5)        ! BDF2-LUSGS
              call BDF2_stage(istep)
            case default
              write(*,*) 'Wrong choice of implicit time stepping method, rk=',rk
              stop
            end select
          endif

        endif
       
      enddo
      
      ! if( (.not. IfSteady) .and. ( (.not. IfExplicit) .and. (implicit_time_method==1)) )  call updatevariables_irk
      
      if(.not. IfSteady) then
        ttime= ttime_n + dt
        call compute_p_out
      endif

      if(istep==1 .or. mod(istep,NScreenOutput)==0)then
        call cpu_time(time1) 
        time_all=time1-time0

        write(*,*) 'istep=',istep,'in processor',MyProc
        write(*,*) 'ttime,dt,cputime:',ttime,dt,time_all
        
        ! write(*,*) 'max production_rate_sum=', maxval(production_rate_sum(1:NE))
        ! write(*,*) 'min production_rate_sum=', minval(production_rate_sum(1:NE))

        ! write(*,*) 'max flux_difference_inv', maxval(flux_difference_nosource(1:NE))
        ! write(*,*) 'min flux_difference_inv', minval(flux_difference_nosource(1:NE))

        ! write(*,*) 'max flux_difference=', maxval(flux_difference(1:NE))
        ! write(*,*) 'min flux_difference=', minval(flux_difference(1:NE))

        ! do i=1,NE
        !   sum_of_mass_fraction(i)= abs(1.-sum(PA(nflow+1:Nvar,i)))
        ! enddo

        ! write(*,*) 'max sum error of mass fraction=', maxval(sum_of_mass_fraction(1:NE))

        ! write(*,*) 'max mass fraction=', maxval(PA(nflow+1:Nvar,1:NE))
        ! write(*,*) 'min mass fraction=', minval(PA(nflow+1:Nvar,1:NE))

        call step_residual(istep,time_all,RoRes)
        ! call check_mass_fraction
      endif
      
      call save_start_solution

      if(mod(istep,NWriteMonitor)==0) call write_monitor_point_value
      if(mod(istep,NFileOutput)==0) call output_vtk(istep/NFileOutput)
      if(mod(istep,NFileBackup)==0) call BackupData
     
      if( istop==1 ) exit

    enddo

    !!=======================
    !! Save data
    !!=======================
    call save_common_boundary_info
    call BackupData

    !!=======================
    !! Final output
    !!=======================
    if(istep==nstep_max+1) then
      !! loop finished
      call output_vtk(nstep_max)
      write(*,*) 'Simulation finished on Processor-',MyProc,'at T=', ttime, 'after istep=', nstep_max, 'time steps'
    else
      !! loop quited
      call output_vtk(istep)
      write(*,*) 'Simulation finished on Processor-',MyProc,'at T=', ttime, 'after istep=', istep, 'time steps'
    endif

    !!=======================
    !! MPI finalization
    !!=======================
    ! call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_Finalize(ierr)
    
end program


!!==========================================================
!!  reconstruction
!!==========================================================
subroutine reconstruction
    use mainvar
    implicit none

    call gradient
    if(limit_ID>0) call limiter_Barth

    return
end subroutine

!!==========================================================
!!  flux & sourceterm computation
!!==========================================================
subroutine flux_sourceterm
    use mainvar
    implicit none
    integer:: i

    hWA(:,:)= 0.
    call flux_lax !H_center !Roe !LLF !RF !lax

    ! write(*,*) 'inviscid flux max/min:',maxval(hWA(:,1:NE)),minval(hWA(:,1:NE))

    do i=1,NE
      flux_difference_nosource(i)= (sum(hWA(Nflow+1:Nvar,i)) - hWA(1,i))/vol(i)
    enddo

    if( IfViscous ) call flux_viscous

    ! write(*,*) 'viscous flux max/min:',maxval(hWA(:,1:NE)),minval(hWA(:,1:NE))

    do i=1,NE
      flux_difference(i)= (sum(hWA(Nflow+1:Nvar,i)) - hWA(1,i))/vol(i)
    enddo
    

    call sourceterm



    return
end subroutine

!!==========================================================
!!      exit pressure
!!==========================================================
subroutine compute_p_out
    use mainvar
    implicit none
    ! double precision, parameter:: A_p= 0.1, f_p= 5000.


    ! p_out= p_static
    p_out= p_static*(1. + A_p*sin(2*pi*f_p*ttime))

    return
end subroutine

subroutine  check_mass_fraction
  use mainvar
  use Parallel_Var
  implicit none
  include 'mpif.h'
  integer:: i,ierr
  double precision:: min_mass_fraction,max_mass_fraction,max_err_sum_mf, &
                     min_mass_fraction_all,max_mass_fraction_all,max_err_sum_mf_all

  min_mass_fraction= minval(PA(nflow+1:Nvar,1:NE))
  max_mass_fraction= maxval(PA(nflow+1:Nvar,1:NE))

  max_err_sum_mf=  -1.0e10

  do i=1,NE
    max_err_sum_mf= max(max_err_sum_mf, abs(sum(PA(nflow+1:Nvar,i))-1.0))
  enddo

  call MPI_AllReduce(min_mass_fraction, min_mass_fraction_all, 1, MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(max_mass_fraction, max_mass_fraction_all, 1, MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  call MPI_AllReduce(max_err_sum_mf, max_err_sum_mf_all, 1, MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

  if(myProc==1) then
    write(*,*) 'minimum mass fraction value:', min_mass_fraction_all
    write(*,*) 'maximum mass fraction value:', max_mass_fraction_all
    write(*,*) 'maximum mass fraction sum error:', max_err_sum_mf_all
  endif

  return
end subroutine
