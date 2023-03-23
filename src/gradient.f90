!!====================================
!!   Lest-squares reconstruction
!!====================================
subroutine gradient
    use mainvar
    implicit none
    double precision :: C(maxnei,Nvar)
    integer:: i,j,k,nei

    call update_ghostcell_variables
   
    do i=1,NE

      C=0.
       
      do nei=1, neigh(i)
        j= NeighC(nei,i)
        
        !! average
        do k=1,Nvar
          C(nei,k)= wdis(nei,i)*( PA(k,j)- PA(k,i) )
        enddo
      enddo

      if(Num_ghostcell(i)>0) then
        do nei=1, Num_ghostcell(i)
          j= Ghostcell(nei,i)    
          !! average
          do k=1,Nvar
            C(nei+neigh(i),k)= wdis(nei+neigh(i),i)*( GPA(k,j)- PA(k,i) )
          enddo
        enddo
      endif

      do k=1,Nvar
      gradC(1:MaxOrder,k,i)= matmul(AMLS(1:MaxOrder,1:KFace(i),i),C(1:KFace(i),k))
      enddo

     
    enddo     

    call Connect_Gradient
 
    return
end subroutine

