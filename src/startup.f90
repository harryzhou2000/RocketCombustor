subroutine startup
      use mainvar
      use Parallel_Var
      implicit none
      integer :: ig
      integer, parameter:: iofile=40

      !!=================================
      !! species info
      !!---------------------------------
      call setup_species_info
      
      !!=============================
      !! thermo property
      !!-----------------------------
      call setup_thermo_property

      !!=============================
      !! transport property
      !!-----------------------------
      call setup_transport_property
      
      !!=============================
      !! diffusion velocity
      !!-----------------------------
      call setup_diffusion_velocity

      !!===================================
      !! chemical kinetics
      !!-----------------------------------
      call setup_chemical_kinetics

      !!=============================
      !! oxidizer inlet
      !!-----------------------------
      call setup_oxidizer_info
      
      !!=============================
      !! fuel inlet
      !!-----------------------------
      call setup_fuel_info


      !!======================================================
      !!  set the total number of variables
      !!  take account of flow, turbulence and species
      !!------------------------------------------------------

      call setup_number_of_variables

      !!=======================================
      !!  Read fow solver setting
      !!---------------------------------------
      
      filename= trim("flow/flow_solver_info.in")
      OPEN(iofile,FILE=filename,STATUS='OLD')

      write(*,*) 'reading flow solver setting ...'

      read(iofile,*)
      read(iofile,*) If_turn_on_combustion
      
      read(iofile,*) 
      read(iofile,*) IfRightDomain

      read(iofile,*) 
      read(iofile,*) NPreviousPar

      read(iofile,*) 
      read(iofile,*) Npar

      read(iofile,*)
      read(iofile,*) unit

      read(iofile,*)
      read(iofile,*) IfSteady, IfViscous, IfReadBackup, IfExplicit

      ! read(iofile,*) 
      ! read(iofile,*) Re_in, Ma_in

      ! read(iofile,*) 
      ! read(iofile,*) p_out_ratio

      ! read(iofile,*) 
      ! read(iofile,*) alpha_st

      read(iofile,*) 
      read(iofile,*) implicit_time_method

      read(iofile,*)
      read(iofile,*) cfl,totaltime,nstep_max, save_start_time, save_end_time

      read(iofile,*)
      read(iofile,*) IfFixDt

      read(iofile,*)
      read(iofile,*) InLoop,InnerLoop,relTol

      write(*,*) 'The stopping creterion for steady-state/unsteady-inner-iteration: Res<Res0*', relTol

      read(iofile,*)
      read(iofile,*) NScreenOutput,NFileOutput,NFileBackup,NWriteMonitor

      read(iofile,*) 
      read(iofile,*) rc_detector

      read(iofile,*) 
      read(iofile,*) n_wbap,limit_ID

      read(iofile,*) 
      read(iofile,*) num_global_monitor_points
      allocate(x_monitor_point(num_global_monitor_points));x_monitor_point=0.
      allocate(y_monitor_point(num_global_monitor_points));y_monitor_point=0.
      read(iofile,*) 
      do ig=1,num_global_monitor_points
            read(iofile,*) x_monitor_point(ig),y_monitor_point(ig)
      enddo

      read(iofile,*) 
      read(iofile,*) p_out_ratio,A_p,f_p

      ! read(iofile,*)
      ! read(iofile,*) k_fre,k_sin

      ! read(iofile,*)
      ! read(iofile,*) delta_sin

      close(iofile)

      !!!!----initial parallel enviroment---
      call InitParallel

      !!!!----set default parameters---
      call InitParam
    
    return
end subroutine