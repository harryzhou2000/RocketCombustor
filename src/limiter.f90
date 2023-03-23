
!!=================================================================
!!   Barth Limiter
!!=================================================================
subroutine Limiter_Barth
    use mainvar
    implicit none
    integer:: i,j,k,nei,kf,iface,n1,n2
    double precision:: PA_min(Nvar), PA_max(Nvar), phi(Nvar)
    double precision:: xy(2), dxy(2), vvh(MaxOrder), face_value, phi_face, phi_mass_fraction
    
    ! call ShockDetect
    
!!---------  Barth limiting -------------------
    do i=1, NE
      
      ! if( Smooth(i)<1. ) cycle


      !! Now we look for the max and min of the variables

      PA_min= PA(:,i)
      PA_max= PA(:,i)
      
      do nei=1,  Neigh(i)
        j= NeighC(nei,i)

        do k=1,Nvar
          PA_min(k)= min(PA_min(k),PA(k,j))
          PA_max(k)= max(PA_max(k),PA(k,j))
        enddo

      enddo

      if(Num_ghostcell(i)>0) then
        do nei=1, Num_ghostcell(i)
          j= Ghostcell(nei,i)    
          do k=1,Nvar
            PA_min(k)= min(PA_min(k),GPA(k,j))
            PA_max(k)= max(PA_max(k),GPA(k,j))
          enddo
        enddo
      endif

      !! we need to restrict mass fraction in range [0,1]
      ! do k=Nflow+1,Nvar
      !   PA_min(k)= max(0.,PA_min(k))
      !   PA_max(k)= min(1.,PA_max(k))
      ! enddo

      !! Now we apply the Barth limiting procedure

      phi(:)= 1.

      do k=1,Nvar
        do kf=1,KFace(i)

          !! Less retrictive
          iface = NeighF(kf,i)

          n1= FNode(1,iface)
          n2= FNode(2,iface)
          xy(:)= 0.5*(Coor(:,n1) + Coor(:,n2))

          !! Original Barth Limiter
          ! xy(:)= Coor(:, N(kf,i))

          !! compute values at the middle point
          dxy(:)= xy(:)- CellXY(:,i)
          call fbaseAll(dxy,vvh,i)
          
          face_value= PA(k,i) + sum(gradC(:,k,i)*vvh(:))

          if(face_value>PA(k,i)) then
            phi_face= min(1., (PA_max(k) - PA(k,i))/(face_value - PA(k,i)) )
          elseif(face_value<PA(k,i)) then
            phi_face= min(1., (PA_min(k) - PA(k,i))/(face_value - PA(k,i)) )
          else
            phi_face= 1.
          endif

          phi(k)= min(phi(k),phi_face)

        enddo

        ! if(phi(k)<1) write(*,*) phi(k)
      enddo

      ! phi_mass_fraction= minval(phi(Nflow+1:Nvar))

      ! phi(Nflow+1:Nvar)= phi_mass_fraction

      ! phi_cell(:,i)= phi(:)

      do k=1,Nvar
        gradC(:,k,i)= phi(k)*gradC(:,k,i)
      enddo
      
    enddo

    ! call MPI_STOP

    !! data transformation
    call Connect_Gradient
    
    return
end subroutine

