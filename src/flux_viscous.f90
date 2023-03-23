!!===========================================
!!  viscous flux computation
!!-------------------------------------------
subroutine flux_viscous
    use mainvar
    implicit none
    integer:: i,j,k,ig,NL,NR,n1,n2
    double precision :: xy(2),dxyL(2),dxyR(2),vvh(MaxOrder), &
                        tgd(2,Nvar),Flu_Vis(Nvar),Flu(Nvar),vvhx(MaxOrder,2)

    double precision:: Yl(Nvar-Nflow+1),Yr(Nvar-Nflow+1),YY(Nvar-Nflow+1),&
                       dYdx(Nvar-Nflow+1),dYdy(Nvar-Nflow+1),&
                       Dk(Nvar-Nflow+1),&
                       u_l(Nvar),u_r(Nvar),entha(Nvar-Nflow+1),spx(Nvar-Nflow+1),spy(Nvar-Nflow+1)

    double precision :: dx,dy,sav,sav1,sav2,sav1n,sav2n,&
              rl,rr,ul,ur,vl,vr,pl,pr,tl,tr,uu,vv,pp,tem,&
              Cpl,gamal,Rcpcvl,Cpr,gamar,Rcpcvr,Cp,gama,Rcpcv,&
              dudx,dudy,dvdx,dvdy,dpdx,dpdy,dtdx,dtdy,dtdn,vmulc,akmu,vlac,div,&
              txx,txy,tyy,btx,bty,dxk,qx,qy,uicx,uicy,dlrx,dlry,vn,ac,&
              snx,sny,stx,sty,dutdx,dutdy,dundx,dundy
    ! double precision,parameter:: yita=0.5
    double precision,parameter:: yita= 4./3
      
    do i=1,  NF
      NL= Neighbor(1,i)
      NR= Neighbor(2,i)
      DX=     VECX(1,i)
      DY=     VECX(2,i)
      n1= FNode(1,i)
      n2= FNode(2,i)
      sav1= dx
      sav2= dy
      sav = max(sqrt(dx*dx+dy*dy),small)
      sav1n=sav1/sav
      sav2n=sav2/sav
        
      do ig=1, NG
      
        xy(:)= Coor(:,n1)+ wxx(ig)*(Coor(:,n2)-Coor(:,n1))
      
        !! Left and Right state
        dxyL(:)= xy(:)- CellXY(:,NL)
        call fbaseAll(dxyL,vvh,NL)
      
        do k=1,Nvar
          u_l(k)=PA(k,NL)+ sum(gradC(:,k,NL)*vvh(:))
        enddo
      
        pl= u_l(1)
        ul= u_l(2)
        vl= u_l(3)
        tl= u_l(4)

        Yl(1:Nvar-Nflow)= u_l(Nflow+1:Nvar)
        Yl(Nvar-Nflow+1)= 1. - sum(Yl(1:Nvar-Nflow))
        
        ! if( tl< 300 .or. tl> 5000)then
        !   write(*,*) 'error in flux_viscous, the left cell is:',NL
        !   write(*,*) 'location of the left cell is:',cellxy(1:2,NL)
        !   write(*,*) 'the temperature is:', tl
        !   call MPI_STOP
        !   endif

        call ComputeGasParameter(Yl, tl, Cpl, gamal, Rcpcvl)
        rl= pl/(Rcpcvl*tl)
      
        ! if( rl<0. .or. pl<0. .or. tl<0. .or. minval(Yl)<0. .or. maxval(Yl)>1.)then
        !   pl= PA(1,NL)
        !   ul= PA(2,NL)
        !   vl= PA(3,NL)
        !   tl= PA(4,NL)
        !   Yl= PA(Nflow+1:Nvar,NL)
      
        !   call ComputeGasParameter(Yl, tl, Cpl, gamal, Rcpcvl)
        !   rl= pl/(Rcpcvl*tl)
        ! endif
      
      
        if( NR>0 )then
          dxyR(:)= xy(:)- (CellXY(:,NR)-SideOfs(:,i))
          call fbaseAll(dxyR,vvh,NR)
      
          do k=1,Nvar
            u_r(k)=PA(k,NR)+ sum(gradC(:,k,NR)*vvh(:))
          enddo
      
          pr= u_r(1)
          ur= u_r(2)
          vr= u_r(3)
          tr= u_r(4)

          Yr(1:Nvar-Nflow)= u_r(Nflow+1:Nvar)
          Yr(Nvar-Nflow+1)= 1. - sum(Yr(1:Nvar-Nflow))
        
          ! if( tr< 300 .or. tr> 5000)then
          ! write(*,*) 'error in flux_viscous, the right cell is:',NR
          ! write(*,*) 'location of the left cell is:',cellxy(1:2,NR)
          ! write(*,*) 'the temperature is:', tr
          ! call MPI_STOP
          ! endif

          call ComputeGasParameter(Yr, tr, Cpr, gamar, Rcpcvr)
          rr= pr/(Rcpcvr*tr)

          ! if(rr<0. .or. pr<0. .or. tr<0. .or. minval(Yr)<0. .or. maxval(Yr)>1. )then
          !   pr= PA(1,NR)
          !   ur= PA(2,NR)
          !   vr= PA(3,NR)
          !   tr= PA(4,NR)
          !   Yr= PA(Nflow+1:Nvar,NR)
          !   call ComputeGasParameter(Yr, tr, Cpr, gamar, Rcpcvr)
          !   rr= pr/(Rcpcvr*tr)
          ! endif

        elseif( NR== -1 )then !! viscous solid wall
          rr= rl; pr= pl; tr= tl
          ul= 0.; vl= 0.
          ur= 0.; vr= 0.
          Yr= Yl; Cpr= Cpl; Rcpcvr= Rcpcvl; gamar= gamal

        elseif( NR== -2 )then !! symmetry
          call removedirection(ul,vl,sav1n,sav2n)
          rr= rl; pr= pl; tr=tl
          ur= ul; vr= vl
          Yr= Yl; Cpr= Cpl; Rcpcvr= Rcpcvl; gamar= gamal

        elseif(NR== -3)then !! oxidizer inlet

          call oxidizer_inlet_bc(pl,pr,ur,vr,tr,rr,Yr,Cpr,gamar,Rcpcvr)

          pl= pr; ul= ur; vl= vr; tl= tr; rl= rr 
          Yl= Yr; Cpl= Cpr; gamal= gamar; Rcpcvl= Rcpcvr 

        elseif(NR== -4)then !! fuel inlet

          call fuel_inlet_bc(pl,pr,ur,vr,tr,rr,Yr,Cpr,gamar,Rcpcvr)

          pl= pr; ul= ur; vl= vr; tl= tr; rl= rr 
          Yl= Yr; Cpl= Cpr; gamal= gamar; Rcpcvl= Rcpcvr

        elseif(NR== -5)then !! nozzle outlet
          
          ac= sqrt(gamal*pl/rl)
          vn= ul*sav1n + vl*sav2n

          if(abs(vn)>ac) then  !! supersonic
            rr= rl
            ur= ul
            vr= vl
            pr= pl
            tr= tl
          else                  !! subsonic
            pr= p_out
            ur= ul + sav1n*(pl-pr)/(rl*ac)
            vr= vl + sav2n*(pl-pr)/(rl*ac)
            rr= rl + (pr-pl)/(ac**2)
            tr= pr/(Rcpcvl*rr)

            pl= pr; ul= ur; vl= vr; tl= tr; rl= rr 
          endif

          Yr= Yl; gamar= gamal; Rcpcvr= Rcpcvl; Cpr= Cpl

        elseif( NR== -6 )then !! inviscid solid wall
          call removedirection(ul,vl,sav1n,sav2n)
          rr= rl; pr= pl; tr=tl
          ur= ul; vr= vl
          Yr= Yl; Cpr= Cpl; Rcpcvr= Rcpcvl; gamar= gamal

        else
          rr= rl; pr= pl; tr= tl
          ur= ul; vr= vl
          Yr= Yl; Cpr= Cpl; Rcpcvr= Rcpcvl; gamar= gamal
        endif
      
    !!------------------------------
    !! viscous flux
    !!------------------------------
        pp=  0.5*(pl+pr )
        uu=  0.5*(ul+ur )
        vv=  0.5*(vl+vr )
        tem= 0.5*(tl+tr )
        YY=  0.5*(Yl+Yr )

        ! if( tem< 300 .or. tem> 5000)then
        !   write(*,*) 'error in flux_viscous, the right cell is:',NR
        !   write(*,*) 'location of the left cell is:',cellxy(1:2,NR)
        !   write(*,*) 'the temperature is tem=:', tem
        !   call MPI_STOP
        ! endif

        call ComputeGasParameter(YY,tem,Cp,gama,Rcpcv)
        rr= pp/(Rcpcv*tem)  !!pp/((gama-1)*rr)
      
        !! NL 
        call fbaseGrad_dxdy( dxyL,vvhx,NL )
        do k =1, Nvar
          tgd(1,k)= sum( gradC(:,k,NL)*vvhx(:,1) )
          tgd(2,k)= sum( gradC(:,k,NL)*vvhx(:,2) )
        enddo
      
        !! NR
        if(NR>0)then
          ! dxk=min(vol(NL),vol(NR))/sav
          
          dlrx= CellXY(1,NL) - (CellXY(1,NR)- SideOfs(1,i))
          dlry= CellXY(2,NL) - (CellXY(2,NR)- SideOfs(2,i))
          dxk= abs(dlrx*sav1n+dlry*sav2n)

          call fbaseGrad_dxdy( dxyR,vvhx,NR )
          !! gradient obtained by dGRP formula
          do k =1, Nvar
            tgd(1,k)= 0.5*(tgd(1,k)+ sum( gradC(:,k,NR)*vvhx(:,1) ) ) + &
                      yita*sav1n*(u_r(k)-u_l(k))/dxk
            tgd(2,k)= 0.5*(tgd(2,k)+ sum( gradC(:,k,NR)*vvhx(:,2) ) ) + &
                      yita*sav2n*(u_r(k)-u_l(k))/dxk
          enddo
        endif
        

        dpdx= tgd(1,1)
        dpdy= tgd(2,1)

        dudx= tgd(1,2)
        dudy= tgd(2,2)
        
        dvdx= tgd(1,3)
        dvdy= tgd(2,3)
        
        dtdx= tgd(1,4)
        dtdy= tgd(2,4)

        dYdx(1:Nvar-Nflow)= tgd(1,Nflow+1:Nvar)
        dYdy(1:Nvar-Nflow)= tgd(2,Nflow+1:Nvar)
        
        if( NR== -1 .or. NR== -6 ) then !! viscous & inviscid solid wall
          call removedirection(dtdx,dtdy,sav1n,sav2n)
        endif

        if( NR== -2 ) then !! symmetry
          snx=  sav1n; sny= sav2n
          stx= -sny; sty= snx

          !! pressure gradient correction
          call removedirection(dpdx,dpdy,snx,sny) 
          !! temperature gradient correction
          call removedirection(dtdx,dtdy,snx,sny) 
          !! mass fraction gradient correction
          do k=1,Nvar-Nflow
            call removedirection(dYdx(k),dYdy(k),snx,sny)
          enddo

          !! velocity gradient correction
          dutdx= dudx*stx + dvdx*sty
          dutdy= dudy*stx + dvdy*sty

          dundx= dudx*snx + dvdx*sny
          dundy= dudy*snx + dvdy*sny

          call removedirection(dutdx,dutdy,snx,sny)
          call removedirection(dundx,dundy,stx,sty)

          dudx= dutdx*stx + dundx*snx
          dudy= dutdy*stx + dundy*snx

          dvdx= dutdx*sty + dundx*sny
          dvdy= dutdy*sty + dundy*sny
        endif

        dYdx(Nvar-Nflow+1)= - sum(dYdx(1:Nvar-Nflow))
        dYdy(Nvar-Nflow+1)= - sum(dYdy(1:Nvar-Nflow))

        !!==========================================
        !! compute enthalpy
        !!------------------------------------------
        call compute_enthalpy(tem,entha)

        !!==========================================
        !! compute transport peroperties
        !!------------------------------------------
        call compute_transport_property(tem,YY,vmulc,akmu)

        !!==========================================
        !! compute diffusion velocity
        !!------------------------------------------
        call compute_diffusion_velocity(pp,tem,YY,Dk)

        !! correction velocity
        uicx= sum(Dk(:)*dYdx(:))
        uicy= sum(Dk(:)*dYdy(:))

        vlac =-2./3.*vmulc
        div= dudx+dvdy
        txx= 2.*vmulc*dudx+vlac*div
        tyy= 2.*vmulc*dvdy+vlac*div
        txy= vmulc*(dudy+dvdx)

        qx= -akmu*dtdx - rr*sum((Dk(:)*dYdx(:)-YY(:)*uicx)*entha(:)) !rr*sum(Vk(1,:)*YY(:)*entha(:))
        qy= -akmu*dtdy - rr*sum((Dk(:)*dYdy(:)-YY(:)*uicy)*entha(:)) !rr*sum(Vk(2,:)*YY(:)*entha(:))
        
        btx= uu*txx + vv*txy - qx
        bty= uu*txy + vv*tyy - qy

        spx= rr*(Dk(:)*dYdx(:) - YY(:)*uicx)
        spy= rr*(Dk(:)*dYdy(:) - YY(:)*uicy)
     
        Flu_Vis(1)= 0.
        Flu_Vis(2)=-(sav1*txx+ sav2*txy)
        Flu_Vis(3)=-(sav1*txy+ sav2*tyy)
        Flu_Vis(4)=-(sav1*btx+ sav2*bty)
        Flu_Vis(Nflow+1:Nvar)= -(sav1*spx(1:Nvar-Nflow) + sav2*spy(1:Nvar-Nflow))

        !!-------------------------------------------------

        Flu(:)= Flu_Vis(:)     !! not /Re_in: for non-dimensionalize, because we have diffusion terms

        ! if(isnan(flu(4)) .or. isnan(flu(5)) .or.isnan(flu(6)) .or. isnan(flu(7)) .or. isnan(flu(8))) then
        ! write(*,*) 'there is nan in viscous flux computation'
        ! write(*,*) 'the flux is:',flu
        ! write(*,*) 'NL=',NL,'NR=',NR
        ! write(*,*) 'qx=',qx,'qy=',qy
        ! write(*,*) 'dtdx=',dtdx,'dtdy=',dtdy
        ! write(*,*) 'akmu=',akmu
        ! write(*,*) 'Dk=',Dk
        ! write(*,*) 'dYdx=',dYdx
        ! write(*,*) 'dYdy=',dYdy
        ! write(*,*) 'YY=',YY
        ! write(*,*) 'uicx=',uicx,'uicy=',uicy
        ! write(*,*) 'entha=',entha
        ! stop
        ! endif

        do k=1, Nvar
        hWA(k,NL)= hWA(k,NL) - Flu(k)*wgg(ig)
        if( NR>0 .and. Fproperty(i)/=100)&
        hWA(k,NR)= hWA(k,NR) + Flu(k)*wgg(ig)
        enddo

      enddo
    
    enddo
      
      return
end subroutine

subroutine removedirection(vx,vy,nx,ny)
    implicit none
    double precision:: vx,vy,vn
    double precision, intent(in):: nx,ny

    vn= vx*nx + vy*ny
    vx= vx-vn*nx
    vy= vy-vn*ny

end subroutine

subroutine reversedirection(vx,vy,nx,ny)
    implicit none
    double precision:: vx,vy,vn
    double precision, intent(in):: nx,ny

    vn= vx*nx + vy*ny
    vx= vx-2*vn*nx
    vy= vy-2*vn*ny

end subroutine
