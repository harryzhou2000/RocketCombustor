!!------------------------------------------------------------------
!! implicit Runge-Kutta method, using dual time stepping approach
!!------------------------------------------------------------------
subroutine trapezoidal_stage(istep)
    use mainvar
    use Parallel_Var
    implicit none
    integer:: inner_rk,istep
    double precision:: RoRes, RoRes_tmp

    do inner_rk=1,InnerLoop

      ! write(*,*) 'inner step=',inner_rk

      call timestep

      call reconstruction
      ! write(*,*) 'gradient max/min:',maxval(gradC(:,:,1:NE)),minval(gradC(:,:,1:NE))
      call flux_sourceterm
      ! write(*,*) 'flux max/min:',maxval(hWA(:,1:NE)),minval(hWA(:,1:NE))
      ! write(*,*) 'flux_sourceterm computation finished at processor', MyProc
      ! stop
      call trapezoidal_LU_SGS(Inner_rk)
      ! write(*,*) 'cell average max/min:',maxval(PA(:,1:NE)),minval(PA(:,1:NE))
      ! stop


      ! write(*,*) 'the smallest temperature', minval(PA(4,1:NE)), &
      !            ' is in location:', cellxy(:,minloc(PA(4,1:NE))), 'in processor', MyProc


      call inner_iteration_residual(istep,inner_rk,RoRes)
      if( inner_rk==1 ) RoRes_tmp= RoRes

      !! check convergence
      ! if( log10(RoRes_tmp)- log10(RoRes) > OrderOfResidual .or.  inner_rk==InnerLoop) then
      if( RoRes < RoRes_tmp*relTol .or.  inner_rk==InnerLoop) then      
        call reconstruction
        call flux_sourceterm
        ! hWAP(rk,:,:)=hWA(:,:)

        call save_common_boundary
        if(IfRightDomain) call update_right_solution

        exit
      endif

      if( isnan(WA(1,1)) )then
        write(*,*) 'there is NAN in step ',istep,'innerloop',inner_rk,'in processor', MyProc
        call MPI_STOP
      endif

    enddo

    return
end subroutine

    subroutine trapezoidal_LU_SGS(Inner_rk)
    use mainvar
    use species
    implicit none
    integer :: inner_rk
    integer:: i,j,k,iface,nei,NL,NR,it
    double precision :: dflux(Nvar),yy(Nvar-Nflow+1),ry(Nvar-Nflow+1),dyy(Nvar-Nflow+1),dry(Nvar-Nflow+1),&
              entha(Nvar-Nflow+1),hh
    double precision,allocatable:: dd(:)
    double precision :: Res,dx,dy,rr,rou,rov,roe,uu,vv,pp,tt,&
              drr,dru,drv,dre,duu,dvv,dpp,dtt,unormal,dunormal,&
              Cp, gama, Rcpcv
    double precision:: JacoMat(Nvar,Nvar),dwdp(Nvar-Nflow,Nvar),rhs(Nvar),sol(Nvar)

    allocate(dd(NE))

    !! at the start of an inner iteration
    if( inner_rk==1 )then  
      do i=1,NE
        fixed(:,i)=vol(i)/dt*WA(:,i) + 0.5*hWA(:,i)
      enddo
    endif
    
     !! diagonal value             
    dd(:) =0.
    do i=1, NE 
    do nei=1,KFACE(i)  !! Tri/Quad= 3/4
      iface= NeighF(nei,i)
      dd(i)= dd(i)+ 0.5*0.5*alpaF(iface)
    enddo
    dd(i)= dd(i)+ vol(i)*(1./step(i)+ 1./dt)
    enddo
    
    !! update residual
    do i=1, NE
    do k=1, Nvar
      !! Res is right hand side
      Res= 0.5*hWA(k,i)- vol(i)/dt*WA(k,i)+fixed(k,i)
      hWA(k,i)= Res
    enddo
    enddo
    
    !! forward sweep
    do i=1, NE

      ! if(PA(4,i)<0) then
      !       write(*,*) 'Wring temperature, tt=',tt,'of cell i=',i
      !       stop
      !     endif

      !! compute Jacobian Matrix
      call Compute_Con_to_Pri_Jaco(PA(:,i),JacoMat)
      JacoMat= JacoMat*dd(i)

      call compute_derivative_production_rate(PA(:,i),dwdp)

      do j=1,Nvar-Nflow
        JacoMat(nflow+j,:)= JacoMat(nflow+j,:) - 0.5*vol(i)*dwdp(j,:)
      enddo

      !! collect coef for every Delta_U(j)
      dflux(:)= 0.
      do nei=1,KFACE(i)  !! Tri/Quad= 3/4
        iface= NeighF(nei,i)
        dx= VECX(1,iface)   !! (dx,dy): normal vector
        dy= VECX(2,iface)
        NR= Neighbor(2,iface)
        if( i==NR )then   !! assure NL=i, maybe i is not left cell of iface
          NR= Neighbor(1,iface)
          dx= -dx  !!??
          dy= -dy
        endif

        if( NR>0 .and. NR<i )then

          pp= PA(1,NR)
          uu= PA(2,NR)
          vv= PA(3,NR)
          tt= PA(4,NR)
          yy(1:Nvar-Nflow)= PA(nflow+1:Nvar,NR); yy(Nvar-Nflow+1)= 1. - sum(yy(1:Nvar-Nflow))

          ! if(tt<0) then
          !   write(*,*) 'Wring temperature, tt=',tt,'of cell NR=',NR
          !   stop
          ! endif
          
          call ComputeGasParameter(yy, tt, Cp, gama, Rcpcv)
          
          call compute_enthalpy(tt,entha)
          hh= sum(yy(:)*entha(:))

          rr = pp/(Rcpcv*tt)
          rou= rr*uu
          rov= rr*vv
          roe= rr*(hh + 0.5*(uu**2+vv**2) ) - pp
          ry = rr*yy

          dpp= hWA(1,NR)
          duu= hWA(2,NR)
          dvv= hWA(3,NR)
          dtt= hWA(4,NR)
          dyy(1:Nvar-Nflow)= hWA(Nflow+1:Nvar,NR); dyy(Nvar-Nflow+1)= -sum(dyy(1:Nvar-Nflow))

          drr= (dpp - rr*Rcpcv*dtt - rr*tt*sum(R_species(:)*dyy(:)))/(Rcpcv*tt)
          dru= rr*duu + uu*drr
          drv= rr*dvv + vv*drr
          dre= rr*Cp*dtt + rr*sum(entha(:)*dyy(:)) + hh*drr + &
                  rr*(uu*duu + vv*dvv) + 0.5*(uu**2+vv**2)*drr - dpp
          dry(:)= rr*dyy(:) + yy(:)*drr
          
          unormal =  uu*dx+  vv*dy
          dunormal= duu*dx+ dvv*dy
          
          ! ( F(Uj+dUj)-F(Uj) )*Sj, only use cell j

          dflux(1)=dflux(1)+0.5*(-alpaF(iface)*drr+ dru*dx+drv*dy )
          dflux(2)=dflux(2)+0.5*(-alpaF(iface)*dru+ dru*unormal+ rou*dunormal+ dpp*dx )
          dflux(3)=dflux(3)+0.5*(-alpaF(iface)*drv+ drv*unormal+ rov*dunormal+ dpp*dy )
          dflux(4)=dflux(4)+0.5*(-alpaF(iface)*dre+ dunormal*(roe+pp)+ unormal*(dre+dpp) )
          dflux(nflow+1:Nvar)= dflux(nflow+1:Nvar) + 0.5*(-alpaF(iface)*dry(1:Nvar-Nflow) + &
            dry(1:Nvar-Nflow)*unormal + ry(1:Nvar-Nflow)*dunormal)
      
       endif
      enddo

      !! compute d_PA(:,i)
      rhs(:)= hWA(:,i)- 0.5*dflux(:)
      call advanced_linear_solve_lapack(JacoMat,Nvar,rhs,sol)
      hWA(:,i)= sol

      ! do k=1,Nvar
      !   hWA(k,i)= (hWA(k,i)- asj(rk,rk)*dflux(k))/dd(i)
      ! enddo
    enddo

    ! stop
    
    
    !!!------- data transform  -------

       call Connect_LUSGS
       
       
    !! reverse sweep
    do i=NE, 1, -1
      
      !! compute Jacobian Matrix
      call Compute_Con_to_Pri_Jaco(PA(:,i),JacoMat)
      JacoMat= JacoMat*dd(i)

      call compute_derivative_production_rate(PA(:,i),dwdp)

      do j=1,Nvar-Nflow
        JacoMat(nflow+j,:)= JacoMat(nflow+j,:) - 0.5*vol(i)*dwdp(j,:)
      enddo

      !! collect coef for every Delta_U(j)
      dflux(:)= 0.
      do nei=1, KFACE(i)  !! Tri/Quad= 3/4
        iface= NeighF(nei,i)
        DX= VECX(1,iface)   !! (dx,dy): normal vector
        DY= VECX(2,iface)
        NR= Neighbor(2,iface)
        if( i==NR )then   !! assure NL=i, maybe i is not left cell of iface
          NR= Neighbor(1,iface)
          dx= -dx  !!??
          dy= -dy
        endif
        
        if( NR>0 .and. NR>i )then
        
          pp= PA(1,NR)
          uu= PA(2,NR)
          vv= PA(3,NR)
          tt= PA(4,NR)
          yy(1:Nvar-Nflow)= PA(nflow+1:Nvar,NR); yy(Nvar-Nflow+1)= 1. - sum(yy(1:Nvar-Nflow))
          
          call ComputeGasParameter(yy, tt, Cp, gama, Rcpcv)
          
          call compute_enthalpy(tt,entha)
          hh= sum(yy(:)*entha(:))

          rr = pp/(Rcpcv*tt)
          rou= rr*uu
          rov= rr*vv
          roe= rr*(hh + 0.5*(uu**2+vv**2) ) - pp
          ry = rr*yy

          dpp= hWA(1,NR)
          duu= hWA(2,NR)
          dvv= hWA(3,NR)
          dtt= hWA(4,NR)
          dyy(1:Nvar-Nflow)= hWA(Nflow+1:Nvar,NR); dyy(Nvar-Nflow+1)= -sum(dyy(1:Nvar-Nflow))

          drr= (dpp - rr*Rcpcv*dtt - rr*tt*sum(R_species(:)*dyy(:)))/(Rcpcv*tt)
          dru= rr*duu + uu*drr
          drv= rr*dvv + vv*drr
          dre= rr*Cp*dtt + rr*sum(entha(:)*dyy(:)) + hh*drr + &
                  rr*(uu*duu + vv*dvv) + 0.5*(uu**2+vv**2)*drr - dpp
          dry(:)= rr*dyy(:) + yy(:)*drr
          
          unormal =  uu*dx+  vv*dy
          dunormal= duu*dx+ dvv*dy
          
          ! ( F(Uj+dUj)-F(Uj) )*Sj, only use cell j

          dflux(1)=dflux(1)+0.5*(-alpaF(iface)*drr+ dru*dx+drv*dy )
          dflux(2)=dflux(2)+0.5*(-alpaF(iface)*dru+ dru*unormal+ rou*dunormal+ dpp*dx )
          dflux(3)=dflux(3)+0.5*(-alpaF(iface)*drv+ drv*unormal+ rov*dunormal+ dpp*dy )
          dflux(4)=dflux(4)+0.5*(-alpaF(iface)*dre+ dunormal*(roe+pp)+ unormal*(dre+dpp) )
          dflux(nflow+1:Nvar)= dflux(nflow+1:Nvar) + 0.5*(-alpaF(iface)*dry(1:Nvar-Nflow) + &
            dry(1:Nvar-Nflow)*unormal + ry(1:Nvar-Nflow)*dunormal)
        
       endif
      enddo

      !! compute d_PA(:,i)
      rhs(:)= 0.5*dflux(:)
      call advanced_linear_solve_lapack(JacoMat,Nvar,rhs,sol)
      hWA(:,i)= hWA(:,i) - sol
      
      !! compute d_WA(4,i)
      ! do k=1,Nvar
      !   hWA(k,i)= hWA(k,i)- asj(rk,rk)*dflux(k)/dd(i)
      ! enddo
    enddo
    
    ! write(*,*) 'p max, min=', maxval(PA(1,:NE)),minval(PA(1,:NE))
    ! write(*,*) 'u max, min=', maxval(PA(2,:NE)),minval(PA(2,:NE))
    ! write(*,*) 'v max, min=', maxval(PA(3,:NE)),minval(PA(3,:NE))
    ! write(*,*) 'T max, min=', maxval(PA(4,:NE)),minval(PA(4,:NE))
    ! write(*,*) 'Y_CH4 max, min=', maxval(PA(5,:NE)),minval(PA(5,:NE))
    ! write(*,*) 'Y_O2 max, min=', maxval(PA(6,:NE)),minval(PA(6,:NE))
    ! write(*,*) 'Y_CO2 max, min=', maxval(PA(7,:NE)),minval(PA(7,:NE))
    
    !! update primitive variables
    do i=1,NE

      PA(:,i)= PA(:,i)+ hWA(:,i)
      ! if(PA(4,i)<0) then
      !   write(*,*) 'PA=',PA(:,i)
      !   write(*,*) 'hWA=',hWA(:,i)
      !   write(*,*) 'cell center=',cellxy(:,i)
      ! endif
    enddo

    ! stop

    call Connect_CellAverage

    ! stop

    call primitve_to_conserved_variables

    ! do i=1,NE
    ! do k=1,Nvar
    !   WA(k,i)= WA(k,i)+ hWA(k,i)
    ! enddo
    ! enddo
    
    ! call update_primitve_variables
  

    deallocate(dd)

  return
end subroutine
