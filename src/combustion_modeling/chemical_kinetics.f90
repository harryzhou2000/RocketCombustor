!!===============================
!! Chemical Kinetics
!!-------------------------------
module chemical_kinetics
    implicit none
    integer:: code_chemical_model
    integer:: Nreact
    double precision,allocatable:: coef_f_react(:,:),coef_r_react(:,:),&
                                   Af_react(:),beta_react(:),E_react(:),&
                                   dentropy_react(:), denthalpy_react(:)
    double precision:: p_a

end module


subroutine setup_chemical_kinetics
    use species
    use chemical_kinetics
    implicit none
    integer:: i,j
    character*50:: filename
    integer,parameter :: iofile=111

    filename= trim("combustion/chemical_kinetics.in")
    OPEN(iofile,FILE=filename,STATUS='OLD')

    write(*,*) 'reading chemical kinetics info ...'
    
    read(iofile,*) 
    read(iofile,*) Nreact, code_chemical_model

    call allocate_chemical_kinetics_array

    read(iofile,*) ! reaction Konstant coefficients f/r
    do i=1,Nreact
        do j=1,Nsp
            read(iofile,*) coef_f_react(j,i),coef_r_react(j,i)
        enddo
    enddo

    read(iofile,*) ! reaction k's formula constants Af beta E
    do i=1,Nreact
        read(iofile,*) Af_react(i),beta_react(i),E_react(i)
    enddo

    read(iofile,*) ! ? used ?
    do i=1,Nreact
        read(iofile,*) dentropy_react(i), denthalpy_react(i)
    enddo

    read(iofile,*)
    read(iofile,*) p_a
   
    close(iofile)

    return
end subroutine

subroutine allocate_chemical_kinetics_array
    use species
    use chemical_kinetics
    implicit none

    allocate(coef_f_react(Nsp,Nreact));coef_f_react=0.
    allocate(coef_r_react(Nsp,Nreact));coef_r_react=0.

    allocate(Af_react(Nreact));Af_react=0.
    allocate(beta_react(Nreact));beta_react=0.
    allocate(E_react(Nreact));E_react=0.

    allocate(dentropy_react(Nreact));dentropy_react=0.
    allocate(denthalpy_react(Nreact));denthalpy_react=0.

    return
end subroutine


! subroutine compute_production_rate(T,rho,Y,production_rate)
!     use species
!     use chemical_kinetics
!     implicit none
!     double precision:: T,rho,Y(Nsp),react_rate(Nreact),production_rate(Nsp)
!     double precision:: kf_react,kc_react,kr_react
!     integer::j,ks

!    	do j=1,Nreact
!         kf_react= Af_react(j)*T**beta_react(j)*exp(-E_react(j)/(Rconstant*T))

!         kc_react= (p_a/(Rconstant*T))**sum(coef_r_react(:,j)-coef_f_react(:,j))&
!                   *exp(dentropy_react(j)/Rconstant - denthalpy_react(j)/(Rconstant*T))

!         kr_react= kf_react/kc_react

!         do ks=1,Nsp
!           kf_react= kf_react*(Y(ks)*rho/mole_weight(ks))**coef_f_react(ks,j)
!           kr_react= kr_react*(Y(ks)*rho/mole_weight(ks))**coef_r_react(ks,j)
!         enddo

!         react_rate(j)= kf_react - kr_react
!      enddo

!      do ks=1,Nsp
!         production_rate(ks)= mole_weight(ks)*sum((coef_r_react(ks,:)-coef_f_react(ks,:))*react_rate(:))
!      enddo

!     return
! end subroutine

subroutine compute_production_rate(T,rho,Y2,production_rate) ! ! warning, implemented specially
    use species
    use chemical_kinetics
    implicit none
    double precision:: T,rho,Y(Nsp),Y2(Nsp),react_rate(Nreact),production_rate(Nsp)
    double precision:: kf_react,kc_react,kr_react
    integer::j,ks

    Y= Y2
    Y= max(0.,Y)
    Y= min(1.,Y)

    

    if (code_chemical_model == 0) then
        ! do j=1,Nreact

        kf_react= Af_react(1)*exp(-E_react(1)/(Rconstant*T))

        !write(*,*) 'kf_react=',kf_react

        kf_react= kf_react*(Y(1)*rho/(mole_weight(1)*1000))**0.2
        !write(*,*) 'kf_react=',kf_react
        kf_react= kf_react*(Y(2)*rho/(mole_weight(2)*1000))**1.3
        !write(*,*) 'kf_react=',kf_react

        react_rate(1)= kf_react

        !write(*,*) Af_react(j),T,beta_react(j),E_react(j),Rconstant
        !write(*,*) Y(1),Y(2),rho,mole_weight(1),mole_weight(2)
        ! write(*,*) 'kf_react=' kf_react
    ! enddo
    elseif (code_chemical_model == 1) then
        kf_react = 0.0
        if (T > 710.0) then
            kf_react = 1e4
            ! write(*,*) "T hit"
        endif
        kf_react = kf_react * (Y(1)*rho/(mole_weight(1)*1000)) ** 1.0
        kf_react = kf_react * (Y(2)*rho/(mole_weight(2)*1000)) ** 0.0
        ! write(*,*) "Y1 Y2 mole den"
        ! write(*,*) Y(1)*rho/(mole_weight(1)*1000)
        ! write(*,*) Y(2)*rho/(mole_weight(2)*1000)
        react_rate(1)= kf_react
    else
        write(*,*) "No such code_chemical_model"

        stop 
    endif

     do ks=1,Nsp
        production_rate(ks)= 1000*mole_weight(ks)*(coef_r_react(ks,1)-coef_f_react(ks,1))*react_rate(1)
     enddo
    !  write(*,*) kf_react, production_rate

    return
end subroutine

subroutine compute_derivative_production_rate(pri,dwdp)
    use species
    use mainvar
    use chemical_kinetics
    implicit none
    double precision:: pri(Nvar),dwdp(Nsp-1,Nvar)
    double precision:: p,u,v,T,rho,Y(Nsp-1),Y2(Nsp),Cp,gama,Rcpcv,&
                       af,r1,r2,dr1dr,dr1dy1,dr2dr,dr2dy2,drdp,drdT,&
                       drdy(Nsp-1),dkfdr,dkfdp,dkfdu,dkfdv,dkfdT,dkfdy(Nsp-1),&
                       dkf(Nvar),dter1_dr
    integer::i
    double precision,parameter:: epsilon=1.0e-12

    p= pri(1)
    u= pri(2)
    v= pri(3)
    T= pri(4)
    y= pri(Nflow+1:Nvar)

    y2(1:Nsp-1)= y; y2(Nsp)= 1. - sum(y)
    call ComputeGasParameter(y2, T, Cp, gama, Rcpcv)
    rho= p/(Rcpcv*T)

    Y= max(0.,Y)
    Y= min(1.,Y)

    ! Y(1)= max(Y(1),epsilon)
    if(Y(1)<epsilon) then
        dkf=  0.
    else

        af= Af_react(1)*exp(-E_react(1)/(Rconstant*T))
        r1= Y(1)*rho/(mole_weight(1)*1000)
        dr1dr=  Y(1)/(mole_weight(1)*1000)
        dr1dy1=  rho/(mole_weight(1)*1000)
        
        r2= Y(2)*rho/(mole_weight(2)*1000)
        dr2dr=  Y(2)/(mole_weight(2)*1000)
        dr2dy2=  rho/(mole_weight(2)*1000)

        drdp= rho/p
        drdT= -rho/T
        do i=1,Nsp-1
        drdy(i)= -(R_species(i)-R_species(Nsp))*rho/Rcpcv
        enddo

        if(code_chemical_model == 0) then

            dter1_dr=  0.2*rho**(-0.8)*dr1dr**(0.2)  !0.2*(r1**(-0.8))*dr1dr

            dkfdr= af*( (r1**0.2)*1.3*(r2**0.3)*dr2dr +  dter1_dr*(r2**1.3) )

            dkfdp= dkfdr*drdp
            dkfdu= 0.
            dkfdv= 0.
            dkfdT= af*(r1**0.2)*(r2**1.3)*E_react(1)/(Rconstant*T*T) + dkfdr*drdT
            
            dkfdy(1)=  af*0.2*(r1**(-0.8))*dr1dy1*(r2**1.3) + dkfdr*drdy(1)
            dkfdy(2)=  af*(r1**0.2)*1.3*(r2**0.3)*dr2dy2 + dkfdr*drdy(2)
            dkfdy(3:)= dkfdr*drdy(3:)

        elseif(code_chemical_model == 1) then
            af = 0.0
            if (T > 710.0) then
                af = 1e4
            endif

            dter1_dr=  1*rho**(0)*dr1dr**(1)

            dkfdr= af*( (r1**1)*0*(r2**0)*dr2dr +  dter1_dr*(r2**0) )

            dkfdp= dkfdr*drdp
            dkfdu= 0.
            dkfdv= 0.
            dkfdT= 0.
            
            dkfdy(1)=  af*1*(r1**(0))*dr1dy1*(r2**0) + dkfdr*drdy(1)
            dkfdy(2)=  af*(r1**1)*0*(r2**0)*dr2dy2 + dkfdr*drdy(2)
        else
            write(*,*) "No such code_chemical_model"
            stop 
        endif

        dkf=(/ dkfdp,dkfdu,dkfdv,dkfdT,dkfdy/)
    endif

    ! if(isnan(sum(dkf))) then
    !     write(*,*) 'dkf=',dkf
    !     write(*,*) 'pri=',pri
    !     write(*,*) 'rho=',rho
    !     write(*,*) 'Cp,gama,Rcpcv=',Cp,gama,Rcpcv
    !     write(*,*) 'r1,r2,dkfdr,af',r1,r2,dkfdr,af
    ! endif
    do i=1,Nsp-1
        dwdp(i,:)=1000*mole_weight(i)*(coef_r_react(i,1)-coef_f_react(i,1))*dkf
    enddo

    return
end subroutine

! subroutine compute_dprod(pa,dpa,dw)
!     use species
!     use mainvar
!     use chemical_kinetics
!     implicit none
!     double precision:: pa(Nvar),dpa(Nvar),dw(Nsp-1)
!     double precision:: p,u,v,T,rho,Y(Nsp-1),Y2(Nsp),Cp,gama,Rcpcv,&
!                        af,r1,r2,dkf,dp,du,dv,dT,dy(Nsp),drho,daf,dr1,dr2
!     integer::i

!     p= pa(1)
!     u= pa(2)
!     v= pa(3)
!     T= pa(4)
!     y= pa(Nflow+1:Nvar)

!     y2(1:Nsp-1)= y; y2(Nsp)= 1. - sum(y)
!     call ComputeGasParameter(y2, T, Cp, gama, Rcpcv)
!     rho= p/(Rcpcv*T)

!     Y= max(0.,Y)
!     Y= min(1.,Y)

!     dp= dpa(1)
!     du= dpa(2)
!     dv= dpa(3)
!     dT= dpa(4)
!     dy(1:Nsp-1)= dpa(Nflow+1:); dy(Nsp)= -sum(dy(1:Nsp-1))

!     af= Af_react(1)*exp(-E_react(1)/(Rconstant*T))
!     r1= Y(1)*rho/(mole_weight(1)*1000)
!     r2= Y(2)*rho/(mole_weight(2)*1000)

!     drho= (dp - rho*Rcpcv*dT - rho*T*sum(R_species(:)*dy(:)))/(Rcpcv*T)

!     daf= af*E_react(1)/(Rconstant*T*T)*dT
!     dr1= (Y(1)*drho + dy(1)*rho)/(mole_weight(1)*1000)
!     dr2= (Y(2)*drho + dy(2)*rho)/(mole_weight(2)*1000)

!     dkf= daf*r1**0.2*r2**1.3 + af*0.2*r1**(-0.8)*dr1*r2**1.3 + af*r1**0.2*1.3*r2**0.3*dr2

    
!     do i=1,Nsp-1
!         dw(i)=1000*mole_weight(i)*(coef_r_react(i,1)-coef_f_react(i,1))*dkf
!     enddo

!     return
! end subroutine
