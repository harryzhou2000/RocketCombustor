!!--------------------------------------------------------
!! implicit dual time stepping, with LU-SGS smoother
!!----------------------------------------------------
subroutine lu_sgs_stage(istep)
    use mainvar
    use Parallel_Var
    implicit none
    integer:: istep

      call timestep

      call reconstruction
      ! write(*,*) 'gradient max/min:',maxval(gradC(:,:,1:NE)),minval(gradC(:,:,1:NE))
      call flux_sourceterm
      ! write(*,*) 'flux max/min:',maxval(hWA(:,1:NE)),minval(hWA(:,1:NE))
      ! write(*,*) 'flux_sourceterm computation finished at processor', MyProc
      call LU_SGS_steady
      ! write(*,*) 'cell average max/min:',maxval(PA(:,1:NE)),minval(PA(:,1:NE))


      ! write(*,*) 'the smallest temperature', minval(PA(4,1:NE)), &
      !            ' is in location:', cellxy(:,minloc(PA(4,1:NE))), 'in processor', MyProc

      if( isnan(WA(1,1)) )then
        write(*,*) 'there is NAN in step ',istep,'in processor', MyProc
        call MPI_STOP
      endif

    return
end subroutine

subroutine LU_SGS_steady
    use mainvar
    use species
    implicit none
    integer:: i,j,k,iface,nei,NL,NR
    real*8 :: dflux(Nvar),yy(Nvar-Nflow+1),ry(Nvar-Nflow+1),dry(Nvar-Nflow),entha(Nvar-Nflow+1),hh
    real*8,allocatable:: dd(:)
    real*8 :: dx,dy,rr,rou,rov,roe,uu,vv,pp,tt,&
              drr,dru,drv,dre,duu,dvv,dpp,unormal,dunormal,&
              Cp, gama, Rcpcv
  
    allocate(dd(NE))

    !! diagonal value
    dd(:) =0.
    do i=1, NE
    do nei=1,KFACE(i)  !! Tri/Quad= 3/4
      iface= NeighF(nei,i)
      dd(i)= dd(i)+ 0.5*alpaF(iface)
    enddo
    dd(i)= dd(i)+ vol(i)/step(i)
    enddo
    
    !! forward sweep
    do i=1, NE
      !! collect coef for every Delta_U(j)
      dflux(:)= 0.
      do nei=1,KFACE(i)  !! Tri/Quad= 3/4
        iface= NeighF(nei,i)
        dx= VECX(1,iface)   !! (dx,dy): normal vector
        dy= VECX(2,iface)
        NR= Neighbor(2,iface)
        if( i==NR )then   !! assure NL=i, maybe i is not left cell of iface
          NR= Neighbor(1,iface)
          dx= -dx 
          dy= -dy
        endif

       if( NR>0 .and. NR<i )then
        
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
        
        drr= hWA(1,NR)
        dru= hWA(2,NR)
        drv= hWA(3,NR)
        dre= hWA(4,NR)
        dry(1:Nvar-Nflow)= hWA(nflow+1:Nvar,NR); dry(Nvar-Nflow+1)= drr -sum(dry(1:Nvar-Nflow))

        duu= (dru- uu*drr)/rr
        dvv= (drv- vv*drr)/rr
        dpp= (gama-1)*( dre - sum(entha(:)*dry(:)) + &
                        tt*sum(R_species(:)*dry(:))*gama/(gama-1) &
                        - 0.5*( dru*uu + duu*rou + drv*vv + dvv*rov ))
        
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

      !! compute d_WA(4,i)
      do k=1,Nvar
        hWA(k,i)= (hWA(k,i)- dflux(k))/dd(i)
      enddo
    enddo
    
    !!!------- data transform  -------

       call Connect_LUSGS

    !! reverse sweep
    do i=NE, 1, -1
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
          
          drr= hWA(1,NR)
          dru= hWA(2,NR)
          drv= hWA(3,NR)
          dre= hWA(4,NR)
          dry(1:Nvar-Nflow)= hWA(nflow+1:Nvar,NR); dry(Nvar-Nflow+1)= drr -sum(dry(1:Nvar-Nflow))

          duu= (dru- uu*drr)/rr
          dvv= (drv- vv*drr)/rr
          dpp= (gama-1)*( dre - sum(entha(:)*dry(:)) + &
                          tt*sum(R_species(:)*dry(:))*gama/(gama-1) &
                          - 0.5*( dru*uu + duu*rou + drv*vv + dvv*rov ))
          
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
      
      !! compute d_WA(4,i)
      do k=1,Nvar
        hWA(k,i)= hWA(k,i)- dflux(k)/dd(i)
      enddo
    enddo
    
    !! update physical value
    do i=1,NE
    do k=1,Nvar
      WA(k,i)= WA(k,i)+ hWA(k,i)
    enddo
    enddo

    call update_primitve_variables

    deallocate(dd)
    
    return
end subroutine
