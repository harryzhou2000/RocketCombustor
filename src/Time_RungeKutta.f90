
!!==============================================================
!!              evolution  RK3
!!==============================================================
subroutine erk_stage(istep,rk)
    use mainvar
    use Parallel_Var
    implicit none
    integer:: rk,istep

    if( .not. IfExplicit) call timestep

    call reconstruction
    ! write(*,*) 'gradient max/min:',maxval(gradC(:,:,1:NE)),minval(gradC(:,:,1:NE))
    call flux_sourceterm
    ! write(*,*) 'flux max/min:',maxval(hWA(:,1:NE)),minval(hWA(:,1:NE))
    ! write(*,*) 'flux_sourceterm computation finished at processor', MyProc
    call evolution(rk)
    ! write(*,*) 'cell average max/min:',maxval(PA(:,1:NE)),minval(PA(:,1:NE))


    ! write(*,*) 'the smallest temperature', minval(PA(4,1:NE)), &
    !            ' is in location:', cellxy(:,minloc(PA(4,1:NE))), 'in processor', MyProc


    if( isnan(WA(1,1)) )then
        write(*,*) 'there is NAN in step ',istep,'stage',rk,'in processor', MyProc
        call MPI_STOP
    endif

    return
end subroutine

subroutine evolution(rk)
    use mainvar
    implicit none
    integer :: i,k,rk
    real :: coef(3)

    coef(1:3)= (/ 1., 0.25, 2./3 /)
    ! coet(1:3)= (/ 1., -0.5, 0.5 /)

    if(rk==1)then
      WP(:,:)= WA(:,:)  !! for temporary storage
    endif
    
    do i=1,  NE
    do k=1,  Nvar
      WA(k,i) = (1.-coef(rk))*WP(k,i) + coef(rk)* ( WA(k,i) + step(i)/Vol(i)* hWA (k,i) )
    enddo
    enddo

    call update_primitve_variables
    
    return
end subroutine


