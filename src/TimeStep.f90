!!==============================================================
!!      Time Step
!!==============================================================
subroutine TimeStep
    use mainvar
    implicit none
    include 'mpif.h'
    integer:: i, NL, NR, ierr
    real*8 :: dx,dy,rr,uu,vv,pp,tt,Cp,gama,Rcpcv,yy(Nvar-Nflow+1),&
              unormal, vmul, tem, alpa, alpaVis, dt_all,akmu,prl
    real*8,allocatable:: alpaC(:)

    allocate(alpaC(NE_Vir))
    
    alpaC(:)= 0.
    do i=1, NF
      NL= Neighbor(1,i)
      NR= Neighbor(2,i)
      DX=     VECX(1,i)
      DY=     VECX(2,i)

      !! cell interface 
      if( NR>0 )then

        pp= 0.5*(PA(1,NL) + PA(1,NR))
        uu= 0.5*(PA(2,NL) + PA(2,NR))
        vv= 0.5*(PA(3,NL) + PA(3,NR))
        tt= 0.5*(PA(4,NL) + PA(4,NR))
        yy(1:Nvar-Nflow)= 0.5*(PA(nflow+1:Nvar,NL)+ PA(nflow+1:Nvar,NR))
        yy(Nvar-Nflow+1)= 1. - sum(yy(1:Nvar-Nflow))

        ! if( tt< 300 .or. tt> 5000)then
        !   write(*,*) 'error in computing timestep, the left and right cells are:',NL,NR
        !   write(*,*) 'location of the left cell is:',cellxy(1:2,NL)
        !   write(*,*) 'the temperature is:', tt
        !   call MPI_STOP
        ! endif

        call ComputeGasParameter(yy, tt, Cp, gama, Rcpcv)
        rr= pp/(Rcpcv*tt)

        

      !! boundary face 
      else  
        pp= PA(1,NL) 
        uu= PA(2,NL) 
        vv= PA(3,NL) 
        tt= PA(4,NL)
        yy(1:Nvar-Nflow)= PA(nflow+1:Nvar,NL)
        yy(Nvar-Nflow+1)= 1. - sum(yy(1:Nvar-Nflow))

        ! if(  tt< 300 .or. tt> 5000)then
        !   write(*,*) 'error in computing timestep, the boundary cell is:',NL
        !   write(*,*) 'location of the boundary cell is:',cellxy(1:2,NL)
        !   write(*,*) 'the temperature is:', tt
        !   call MPI_STOP
        ! endif

        call ComputeGasParameter(yy, tt, Cp, gama, Rcpcv)
        rr= pp/(Rcpcv*tt)

        
      endif
      
      if( rr<0. .or. pp<0. .or. tt<0.)then
        write(*,*) 'error in computing timestep, the left and right cells are:',NL,NR
        write(*,*) 'location of the left cell is:',cellxy(1:2,NL)
        call MPI_STOP
      endif
      
      unormal= uu*dx+ vv*dy
      if( IfViscous )then  !! viscous

        call compute_transport_property(tt,yy,vmul,akmu)

        prl= vmul*Cp/akmu
        
        !! consider viscous effect, viscous CFL stability 
        alpa    = abs(unormal)+ sqrt( gama* pp/ rr ) * sqrt(dx*dx+dy*dy)
        alpaVis = 2.*(dx*dx+dy*dy)*max(4./3.,gama/prl)*vmul/rr
        alpaF(i)= alpa+ alpaVis/Vol(NL)  !! ??????
                  alpaC(NL)= alpaC(NL)+ alpa+ alpaVis/Vol(NL)
        if(NR>0)  alpaC(NR)= alpaC(NR)+ alpa+ alpaVis/vol(NR)
        
      else  !! inviscid
        alpa= abs(unormal)+ sqrt( gama* pp/ rr ) * sqrt(dx*dx+dy*dy)
        alpaF(i)= alpa
                  alpaC(NL)= alpaC(NL)+ alpa
        if(NR>0)  alpaC(NR)= alpaC(NR)+ alpa
      endif
    enddo
    
    do i=1, NE
      step(i)= CFL*vol(i)/alpaC(i)
    enddo
    
    deallocate(alpaC)

!!------------------------------------------------------
!! unsteady,  global time step, step(:)= min( step(:) )
!!------------------------------------------------------
    If( IfSteady )return
    if( .not. IfExplicit ) return

    if( IfFixDt) then
      dt= totaltime/nstep_max 
      step(:)= dt
      return
    endif
    
    dt=minval(step(1:NE))
  
    call MPI_AllReduce(dt, dt_all, 1, MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
    dt= dt_all
    step(:)= dt

    return
end subroutine
