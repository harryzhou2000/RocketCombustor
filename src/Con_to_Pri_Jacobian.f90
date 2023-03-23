subroutine Compute_Con_to_Pri_Jaco(pri,jaco)
    use mainvar
    use species
    implicit none
    double precision:: pri(Nvar),jaco(Nvar,Nvar)
    double precision:: p,u,v,t,y(Nsp-1),y2(Nsp),r,entha(Nsp),h,dhdy(Nsp-1),Cp,gama,Rcpcv
    integer:: i,j
    double precision:: drdp,drdu,drdv,drdt,drdy(Nsp-1),&
                       drudp,drudu,drudv,drudt,drudy(Nsp-1),&
                       drvdp,drvdu,drvdv,drvdt,drvdy(Nsp-1),&
                       dredp,dredu,dredv,dredt,dredy(Nsp-1),&
                       drydp(Nsp-1),drydu(Nsp-1),drydv(Nsp-1),drydt(Nsp-1),drydy(Nsp-1,Nsp-1)

    p= pri(1)
    u= pri(2)
    v= pri(3)
    t= pri(4)
    y= pri(Nflow+1:Nvar)

    y2(1:Nsp-1)= y; y2(Nsp)= 1. - sum(y)
    call ComputeGasParameter(y2, t, Cp, gama, Rcpcv)
    r= p/(Rcpcv*t)

    call compute_enthalpy(t,entha)

    h= entha(Nsp)
    do i=1,Nsp-1
    h= h + y(i)*(entha(i)-entha(Nsp))
    dhdy(i)= entha(i)-entha(Nsp)
	enddo

    drdp= 1./(Rcpcv*t)
    drdu= 0.
    drdv= 0.
    drdt= -p/(Rcpcv*t*t)

    do i=1,Nsp-1
    drdy(i)= -(R_species(i)-R_species(Nsp))*p/(t*Rcpcv*Rcpcv)
	enddo

    drudp= u*drdp
    drudu= r 
    drudv= 0.
    drudt= u*drdt
    drudy= u*drdy

    drvdp= v*drdp
    drvdu= 0. 
    drvdv= r
    drvdt= v*drdt
    drvdy= v*drdy

    dredp= (h+0.5*(u*u+v*v))*drdp - 1.
    dredu= r*u
    dredv= r*v
    dredt= r*Cp + (h+0.5*(u*u+v*v))*drdt
    dredy= (h+0.5*(u*u+v*v))*drdy + r*dhdy

    drydp= drdp*y
    drydu= 0.
    drydv= 0.
    drydt= drdt*y

    do i=1,Nsp-1
	    do j=1,i-1
        drydy(i,j)= y(i)*drdy(j)
		enddo
    
        drydy(i,i)= r + y(i)*drdy(i)

	    do j=i+1,Nsp-1
        drydy(i,j)= y(i)*drdy(j)
		enddo
    enddo

    jaco(1,:)=(/drdp,drdu,drdv,drdt,drdy/)
    ! write(*,*) jaco(1,:)
    jaco(2,:)=(/drudp,drudu,drudv,drudt,drudy/)
    ! write(*,*) jaco(2,:)
    jaco(3,:)=(/drvdp,drvdu,drvdv,drvdt,drvdy/)
    ! write(*,*) jaco(3,:)
    jaco(4,:)=(/dredp,dredu,dredv,dredt,dredy/)

    do i=1,Nsp-1
        jaco(i+nflow,:)= (/drydp(i), drydu(i), drydv(i), drydt(i), drydy(i,:)/)
    enddo
    
end subroutine