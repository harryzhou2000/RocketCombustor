!!=============================
!! diffusion velocity
!!-----------------------------
module diffusion_velocity
    implicit none
    double precision:: k_B
    double precision,allocatable:: collision_diameter(:),LJ_energy(:)

    double precision, parameter:: omega_d_A=1.06036,omega_d_B=0.15610,omega_d_C=0.19300,&
                                  omega_d_D=0.47635,omega_d_E=1.03587,omega_d_F=1.52996,&
                                  omega_d_G=1.76474,omega_d_H=3.89411

end module

subroutine setup_diffusion_velocity
    use species
    use diffusion_velocity
    implicit none
    integer:: i
    character*50:: filename
    integer,parameter :: iofile=112

    call allocate_diffusion_array

    filename= trim("combustion/diffusion.in")
    OPEN(iofile,FILE=filename,STATUS='OLD')

    write(*,*) 'reading diffusion info ...'

    read(iofile,*) 
    read(iofile,*) k_B !! Boltzmann's constant
      
    read(iofile,*)
    do i=1,Nsp
      read(iofile,*) collision_diameter(i),LJ_energy(i)
    enddo

    close(iofile)

    return
end subroutine

subroutine allocate_diffusion_array
    use species
    use diffusion_velocity
    implicit none

    allocate(collision_diameter(Nsp));collision_diameter=0.
    allocate(LJ_energy(Nsp));LJ_energy=0.

    return
end subroutine


subroutine compute_diffusion_velocity(pp,tem,YY,Dk)
    use species
    use diffusion_velocity
    implicit none
    double precision,parameter:: epsilon=1.0e-16
    integer:: i,j
    double precision:: pp,tem,YY(Nsp),XX(Nsp),Dk(Nsp),&
                       ww,Dij(Nsp),tstar,omega_d,mw,cd,xdsum

    ww= 1./sum(YY(:)/mole_weight(:))
    XX(:)= ww*YY(:)/mole_weight(:)

    do i=1,Nsp

      do j=1,Nsp
        
        mw= 2./(1./mole_weight(i)+1./mole_weight(j))*1000
        cd= 0.5*(collision_diameter(i)+collision_diameter(j))

        tstar= tem/sqrt(LJ_energy(i)*LJ_energy(j))
        omega_d= omega_d_A*tstar**(-omega_d_B) + omega_d_C*exp(-omega_d_D*tstar) + &
                 omega_d_E*exp(-omega_d_F*tstar) + omega_d_G*exp(-omega_d_H*tstar)   

        Dij(j)=0.0266*sqrt(tem**3) /(pp*sqrt(mw)*cd**2*omega_d)

        ! write(*,*) 'j=',j
        ! write(*,*) 'mw=',mw
        ! write(*,*) 'cd=',cd
        ! write(*,*) 'tstar=',tstar
        ! write(*,*) 'omega_d=',omega_d
        ! write(*,*) 'Dij=',Dij(j)
      enddo

      ! Dk(i)= (1.-XX(i))/(sum(XX(:)/Dij(:))-XX(i)/Dij(i) + epsilon)

      xdsum=0.
      do j=1,i-1
        xdsum= xdsum + XX(j)/Dij(j)
      enddo

      do j=i+1,Nsp
        xdsum= xdsum + XX(j)/Dij(j)
      enddo

      Dk(i)= (1.-XX(i))/(xdsum + epsilon)

      ! write(*,*) Dk(i)
      
      ! if(isnan(Dk(i))) then
      !   write(*,*) 'i=',i
      !   write(*,*) 'Dk(i)=',Dk(i)
      !   write(*,*) 'XX=',XX
      !   write(*,*) 'Dij=',Dij
      !   write(*,*) 'upper=',1.-XX(i)
      !   write(*,*) 'lower=',sum(XX(:)/Dij(:))-XX(i)/Dij(i) + epsilon
      !   write(*,*) 'lower components=',XX(:)/Dij(:), epsilon
      !   write(*,*) 'xdsum=',xdsum, xdsum+epsilon
      !   stop
      ! endif

    enddo

    return
end subroutine