
subroutine sourceterm
    use mainvar
    implicit none
    integer:: i,ks
    real*8 :: source(Nvar),Y(Nvar-Nflow+1),prod(Nvar-Nflow+1)
    real*8 :: rho,u,v,p,t,Cp,gama,Rcpcv

    if(.not. if_turn_on_combustion) return

    do i=1,NE

      p= PA(1,i)
      u= PA(2,i)
      v= PA(3,i)
      t= PA(4,i)
      Y(1:Nvar-Nflow)= PA(Nflow+1:Nvar,i); Y(Nvar-Nflow+1)= 1. - sum(Y(1:Nvar-Nflow))
      !! limit the value of mass fractions to [0,1]
       ! do ks=1,Nvar - Nflow
       !   Y(ks)= max(0., Y(ks)); Y(ks)=min(1.,Y(ks))
       ! enddo

       ! limit mass fraction
      ! Y= limited_mass_fraction(Y)

       ! if( t< 300 .or. t> 5000)then
       !    write(*,*) 'error in sourceterm, the cell is:',i
       !    write(*,*) 'location of the cell is:',cellxy(1:2,i)
       !    write(*,*) 'the temperature is tem=:', t
       !    call MPI_STOP
       !  endif
      
      call ComputeGasParameter(Y, t, Cp, gama, Rcpcv)

      rho= p/(t*Rcpcv)

      !!===================================
      !! Source term computation
      !!-----------------------------------

      source(1:Nflow)= 0.
      
      call compute_production_rate(T,rho,Y,prod)

      source(Nflow+1:Nvar)= prod(1:Nvar-Nflow)

      production_rate(:,i)= prod(1:Nvar-Nflow)

      production_rate_sum(i)= sum(prod(1:Nvar-Nflow))

      hWA(:,i)= hWA(:,i) + vol(i)*source(:)

      ! flux_difference(i)= (sum(hWA(Nflow+1:Nvar,i)) - hWA(1,i))/vol(i)
        
      !write(*,*) source(:) 

       !stop

       ! if(isnan(sum(prod))) then
       !  write(*,*) 'there is nan in sourceterm computation'
       !  write(*,*) 'the sourceterm is:', prod
       !  ! write(*,*) 'NL=',NL,'NR=',NR
       !  ! write(*,*) 'the variables:',rl,rr,ul,ur,vl,vr,pl,pr,acl,acr,coex
       !  stop
       !  endif

    enddo

    return
end subroutine
