subroutine InitFlow_Nozzle
    use mainvar
    implicit none
    integer:: i

    ro_in   = 1.0
    temp_in  = 1.0
    
    Rcpcv  = 1./gama/ (Ma_in**2)
    Cp     = 1./(gama-1.) / (Ma_in**2)
    p_in   = Rcpcv
    
    vx_in  = sqrt(gama*p_in/ro_in)* Ma_in
    vy_in  = 0.

    umu_give= 1./Re_in
    write(*,*) 'umu_give=', umu_give

    write(*,*) 'ro_in,vx_in,vy_in,p_in:',ro_in,vx_in,vy_in,p_in

    p_out= 0.5*p_in  !!0.5,0.75

    write(*,*) 'the static exit pressure is p_out=',p_out
  
    WA(:,:)= 0.
    
    ! do i=1,NE
    !   WA(1,i)= ro_in
    !   WA(2,i)= ro_in* vx_in
    !   WA(3,i)= ro_in* vy_in
    !   WA(4,i)= p_in/(gama-1.)+ 0.5*ro_in*(vx_in**2 + vy_in**2)
    !   ! WA(4,i)= p_out/(gama-1.)+ 0.5*ro_in*(vx_in**2 + vy_in**2)
    ! enddo
    
    do i=1,NE
      WA(1,i)= ro_in
      WA(2,i)= 0 !ro_in* vx_in
      WA(3,i)= 0 !ro_in* vy_in
      ! WA(4,i)= p_in/(gama-1.)
      ! WA(4,i)= p_out/(gama-1.)
      WA(4,i)= (p_in+(CellXY(1,i)+0.3)*(p_out-p_in)/2.9)/(gama-1.)
    enddo

    !! init WP1, WP2,  for implicit scheme
    if( .not. IfExplicit )then
    WP (:,:)= WA(:,:)
    endif
    
    return
end subroutine




subroutine InitFlow_isentropicVortex
  use mainvar
  implicit none
  real*8 :: xv1(2),xv2(2),xv3(2),xv4(2), xy(2),uu(NVar)
  integer:: i,iv
  real*8,external:: GaussInteg_Quad
    
    
    ro_in  = 1.0
    p_in   = 1.0
    vx_in  = 1. !! 0.  !
    vy_in  = 1. !! 0.  !
    Temp_in  = 1.0

    Rcpcv=1.
    Cp= Rcpcv*gama/(gama-1)
    
  do i=1,NE
      xv1= Coor(:,N(1,i))
      xv2= Coor(:,N(2,i))
      xv3= Coor(:,N(3,i))
      xv4= Coor(:,N(4,i))
    
    !! average
    do iv=1,NVar
      WA(iv,i)=  GaussInteg_Quad( fun_WA,xv1,xv2,xv3,xv4 )/vol(i)
    enddo
    
    enddo
    
    !! init WP1, WP2,  for implicit scheme
    if( .not. IfExplicit ) then
    WP (:,:)= WA(:,:)
    endif
    
    
contains
    
    real*8 function fun_WA(xy)
    implicit none
    real*8 :: xx,yy,ss,xp,yp,rad,du,dv,vx,vy,dTT,ro,pp,xy(2)
     
      xx= xy(1); yy= xy(2)

      ss = 5.      
      xp = xx- 5.
      yp = yy- 5.
      rad= sqrt(xp*xp+yp*yp)
      du = ss/2./pi*exp(0.5*(1.-rad*rad)) * (-yp)
      dv = ss/2./pi*exp(0.5*(1.-rad*rad)) *   xp
      vx= vx_in+ du
      vy= vy_in+ dv

      dTT= -(gama-1.)*ss*ss/(8.*gama*pi*pi) * exp(1.-rad*rad)
      ro= (Temp_in+dTT)** (1./(gama-1.))
      pp= ro* (Temp_in+dTT)


      uu(:)=(/ ro, ro*vx, ro*vy, pp/(gama-1)+0.5*ro*(vx**2+vy**2) /)

      fun_WA=uu(iv)
    
     return
     end function
    
end subroutine

subroutine InitFlow_DoubleMach
  use mainvar
  implicit none
  real*8:: xx,yy,flag,rr,uu,vv,pp
  integer::k
    
   
   DO K=1,NE
     
     xx= CellXY(1,k)
     yy= CellXY(2,k)
     flag=sqrt(3.)*(xx-1./6)-yy
     if(flag.gt.0) then
        rr= 1.4
        uu= 0.
        vv= 0.
        pp= 1.0
     else
        rr= 8.0
        uu= 7.1447
        vv=-4.125
        pp= 116.5
     endif

     WA(1,K)= rr
     WA(2,K)= rr*uu
     WA(3,K)= rr*vv
     WA(4,K)= pp/(gama-1.)+ 0.5*rr*(uu**2 + vv**2)
   ENDDO
   
   !! init WP1, WP2,  for implicit scheme
    if( .not. IfExplicit )then
    WP (:,:)= WA(:,:)
    endif
   
end subroutine


subroutine initFlow_viscousshocktube
    use mainvar
    implicit none
    integer:: i
    real*8 :: rr,uu,vv,pp

    Rcpcv=1.
    Cp= Rcpcv*gama/(gama-1)
   
    umu_give= 1./Re_in
    write(*,*) 'umu_give=', umu_give
    
    WA(:,:)= 0.
    do i=1,NE

      if(CellXY(1,i)<0.5)then
        rr= 120.
        uu= 0.
        vv= 0.
        pp= 120./gama
      
      else

        rr= 1.2
        uu= 0.
        vv= 0.
        pp= 1.2/gama
      endif

      WA(1,i)= rr
      WA(2,i)= rr*uu
      WA(3,i)= rr*vv
      WA(4,i)= pp/(gama-1.)+ 0.5*rr*(uu**2 + vv**2)
      
    enddo
    
    !! init WP1, WP2,  for implicit scheme
    if( .not. IfExplicit )then
    WP (:,:)= WA(:,:)
    endif
    
    return
end subroutine