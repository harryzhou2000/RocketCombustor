!!===========================================
!!  Riemann solver
!!---------------------------------
subroutine flux_lax
      use mainvar
      use neuralnetworks
      implicit none
      integer:: i,j,k,NL,NR,n1,n2,ig
      real*8 :: Flu(Nvar),xy(2),dxyL(2),dxyR(2),vvh(MaxOrder),Yl(Nvar-Nflow+1),Yr(Nvar-Nflow+1),entha(Nvar-Nflow+1)
      real*8 :: dx,dy,xmid,ymid,sav1,sav2,sav,sav1n,sav2n,vn,coex, &
                rl,ul,vl,pl,tl,Cpl,gamal,Rcpcvl,hl,el,vaml,acl,ulnormal,rulnormal,&
                rr,ur,vr,pr,tr,Cpr,gamar,Rcpcvr,hr,er,vamr,acr,urnormal,rurnormal,ac

      integer:: id_face,id_face_number
      type(common_bound),pointer:: tcommon_bound
      
      do i=1,  NF
        NL= Neighbor(1,i)
        NR= Neighbor(2,i)
        DX=     VECX(1,i)
        DY=     VECX(2,i)
        n1= FNode(1,i)
        n2= FNode(2,i)
        xmid= 0.5*( Coor(1,n1)+ Coor(1,n2) )
        ymid= 0.5*( Coor(2,n1)+ Coor(2,n2) )
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
        
        pl= PA(1,NL) + sum(gradC(:,1,NL)*vvh(:))
        ul= PA(2,NL) + sum(gradC(:,2,NL)*vvh(:))
        vl= PA(3,NL) + sum(gradC(:,3,NL)*vvh(:))
        tl= PA(4,NL) + sum(gradC(:,4,NL)*vvh(:))

        ! if( tl< 300 .or. tl> 5000)then
        !   write(*,*) 'error in flux_lax, the left cell is:',NL
        !   write(*,*) 'location of the left cell is:',cellxy(1:2,NL)
        !   write(*,*) 'the temperature is:', tl
        !   call MPI_STOP
        ! endif
        
        do k=1,Nvar-Nflow
          Yl(k)= PA(k+Nflow,NL) + sum(gradC(:,k+Nflow,NL)*vvh(:))
        enddo

        Yl(Nvar-Nflow+1)= 1. - sum(Yl(1:Nvar-Nflow))

        call ComputeGasParameter(Yl, tl, Cpl, gamal, Rcpcvl)
        rl= pl/(Rcpcvl*tl)

        ! if(isnan(rl)) then
        !   write(*,*) 'pl average=',PA(1,NL)
        !   write(*,*) 'tl average=',PA(4,NL)
        !   write(*,*) 'gradC=',gradC(:,:,NL)
        !   write(*,*) 'vvh=',vvh
        !   write(*,*) 'pl=',pl
        !   write(*,*) 'tl=',tl
        !   write(*,*) 'Yl=',Yl
        !   write(*,*) 'Rcpcvl=',Yl
        !   stop
        ! endif
        

        ! if( rl<0. .or. pl<0. .or. tl<0. .or. minval(Yl)<0. .or. maxval(Yl)>1.)then
        !   pl= PA(1,NL)
        !   ul= PA(2,NL)
        !   vl= PA(3,NL)
        !   tl= PA(4,NL)
        !   Yl= PA(Nflow+1:Nvar,NL)
        !   call ComputeGasParameter(Yl,tl,Cpl,gamal,Rcpcvl)
        !   rl= pl/(Rcpcvl*tl)
        ! endif 
      
      
      if( NR>0 )then
        dxyR(:)= xy(:)- (CellXY(:,NR)- SideOfs(:,i))
        call fbaseAll(dxyR,vvh,NR)
        
        pr= PA(1,NR) + sum(gradC(:,1,NR)*vvh(:))
        ur= PA(2,NR) + sum(gradC(:,2,NR)*vvh(:))
        vr= PA(3,NR) + sum(gradC(:,3,NR)*vvh(:))
        tr= PA(4,NR) + sum(gradC(:,4,NR)*vvh(:))

        ! if( tr< 300 .or. tr> 5000)then
        !   write(*,*) 'error in flux_lax, the right cell is:',NR
        !   write(*,*) 'location of the left cell is:',cellxy(1:2,NR)
        !   write(*,*) 'the temperature is:', tr
        !   call MPI_STOP
        ! endif
        
        do k=1,Nvar-Nflow
          Yr(k)= PA(k+Nflow,NR) + sum(gradC(:,k+Nflow,NR)*vvh(:))
        enddo

        Yr(Nvar-Nflow+1)= 1. - sum(Yr(1:Nvar-Nflow))

        call ComputeGasParameter(Yr,tr,Cpr,gamar,Rcpcvr)
        rr= pr/(Rcpcvr*tr)
  
        
        ! if(rr<0. .or. pr<0. .or. tr<0. .or. minval(Yr)<0. .or. maxval(Yr)>1.)then
        !   pr= PA(1,NR) 
        !   ur= PA(2,NR) 
        !   vr= PA(3,NR) 
        !   tr= PA(4,NR)
        !   Yr= PA(Nflow+1:Nvar,NR) 
        !   call ComputeGasParameter(Yr,tr,Cpr,gamar,Rcpcvr)
        !   rr= pr/(Rcpcvr*tr)
        ! endif      
      
!!============= boundary face ===============================
       
      elseif( NR== -1 )then !! solid wall

        rr= rl; pr= pl; tr=tl
        ur= -ul ; vr= -vl

        Yr= Yl; gamar= gamal; Rcpcvr= Rcpcvl; Cpr= Cpl

      elseif( NR== -2 )then !! symmetry
        
        rr= rl; pr= pl; tr=tl
        vn= ul*sav1n + vl*sav2n
        ur= ul-2*vn*sav1n
        vr= vl-2*vn*sav2n

        Yr= Yl; gamar= gamal; Rcpcvr= Rcpcvl; Cpr= Cpl

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

        !! supersonic outlet
        ! rr= rl
        ! ur= ul
        ! vr= vl
        ! pr= pl
        ! tr= tl

        !! static pressure 2

        ! ac= sqrt(gamal*pl/rl)

        ! pr= p_out
        ! ur= ul + sav1n*(pl-pr)/(rl*ac)
        ! vr= vl + sav2n*(pl-pr)/(rl*ac)
        ! rr= rl + (pr-pl)/(ac**2)
        ! tr= pr/(Rcpcvl*rr)

        ! pr= 2*pr - pl
        ! ur= 2*ur - ul
        ! vr= 2*vr - vl
        ! tr= 2*tr - tl
        ! rr= 2*rr - rl

        Yr= Yl; gamar= gamal; Rcpcvr= Rcpcvl; Cpr= Cpl

      elseif( NR== -6 )then !! inviscid solid wall
        
        rr= rl; pr= pl; tr=tl
        vn= ul*sav1n + vl*sav2n
        ur= ul-2*vn*sav1n
        vr= vl-2*vn*sav2n

        Yr= Yl; gamar= gamal; Rcpcvr= Rcpcvl; Cpr= Cpl

      elseif(NR<-1000) then !! common boundary

        if(index_common_boundary(1,i)==0 .or. index_common_boundary(1,i)==0) then
          write(*,*) 'error in finding the common boundary face'
          call MPI_STOP
        endif

        id_face= index_common_boundary(1,i)
        id_face_number= index_common_boundary(2,i)

        tcommon_bound => Common_bound_targets(id_face) 

        pr= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+1)
        ur= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+2)
        vr= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+3)
        tr= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+4)

        Yr(1:Nvar-Nflow)= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+Nflow+1:Nvar*id_face_number)
        Yr(Nvar-Nflow+1)= 1. - sum(Yr(1:Nvar-Nflow))

        call ComputeGasParameter(Yr,tr,Cpr,gamar,Rcpcvr)
        rr= pr/(Rcpcvr*tr)
 
      else
        write(*,*) 'error,no face of this type!',i,NR
        stop
      endif

!!===========================================================
      
      if( rl<0. .or. pl<0. .or. (NR>0 .and. (rr<0..or. pr<0.)) )then
        write(*,*) 'error,if:', i, FNode(1,i),FNode(2,i), rl,rr,pl,pr
        stop
      endif
   
      call compute_enthalpy(tl,entha)
      hl= sum(Yl(:)*entha(:))

      call compute_enthalpy(tr,entha)
      hr= sum(Yr(:)*entha(:))
      
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=hl+0.5*vaml
      hr=hr+0.5*vamr
      el=hl-pl/rl
      er=hr-pr/rr

      acl=sqrt(gamal*pl/rl)
      acr=sqrt(gamar*pr/rr)
      
      ! spectrum radius 
       
      ulnormal=(sav1*ul+sav2*vl )
      urnormal=(sav1*ur+sav2*vr )       
      
      coex=max(abs(ulnormal)+acl*sav,abs(urnormal) +acr*sav)

      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      
!     flux terms
      
      Flu(1)= 0.5*(rulnormal+rurnormal-coex*(rr-rl))
      Flu(2)= 0.5*(rulnormal*ul+rurnormal*ur+sav1*(pl+pr)-coex*(rr*ur-rl*ul))
      Flu(3)= 0.5*(rulnormal*vl+rurnormal*vr+sav2*(pl+pr)-coex*(rr*vr-rl*vl))
      Flu(4)= 0.5*(rulnormal*hl+rurnormal*hr-coex*(rr*er-rl*el))

      Flu(Nflow+1:Nvar)= 0.5*(rulnormal*Yl(1:Nvar-Nflow) + rurnormal*Yr(1:Nvar-Nflow) &
                            - coex*(rr*Yr(1:Nvar-Nflow)-rl*Yl(1:Nvar-Nflow)))

      ! if(isnan(flu(4)) .or. isnan(flu(5)) .or.isnan(flu(6)) .or. isnan(flu(7)) .or. isnan(flu(8))) then
      !   write(*,*) 'there is nan in inviscid flux computation'
      !   write(*,*) 'the flux is:',flu
      !   write(*,*) 'NL=',NL,'NR=',NR
      !   write(*,*) 'the variables:',rl,rr,ul,ur,vl,vr,pl,pr,acl,acr,coex
      !   stop
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
