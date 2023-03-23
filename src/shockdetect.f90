!!=================================================================
!!   shock  detector
!!=================================================================
subroutine ShockDetect
    use mainvar
    implicit none
    integer:: i,j,nei,k
    real*8 :: romax,sum1,p,u,v,tem,rho,y(Nvar-Nflow+1),Cp,gama,Rcpcv
    real*8 :: ros(0:MaxNei),dxy(2),vvh(MaxOrder)
    
    smooth=0.
    
    !! use density to detect
    do i=1, NE

      romax= 0.

      do nei=0, Neigh(i)  !! 3
        j= NeighC(nei,i)
        dxy(:)=CellXY(:,i)- (CellXY(:,j)- Ofs(:,nei,i))
        call fbaseAll(dxy,vvh,j)    

        p=   PA(1,j) + sum(gradC(:,1,j)*vvh(:))
        u=   PA(2,j) + sum(gradC(:,2,j)*vvh(:))
        v=   PA(3,j) + sum(gradC(:,3,j)*vvh(:))
        tem= PA(4,j) + sum(gradC(:,4,j)*vvh(:))
        do k=1,Nvar-Nflow
          y(k)= PA(k+nflow,j) + sum(gradC(:,k+nflow,j)*vvh(:))
        enddo
        y(Nvar-Nflow+1)= 1. - sum(y(1:Nvar-Nflow))
        call ComputeGasParameter(y, tem, Cp, gama, Rcpcv)
        rho= p/(Rcpcv*tem)

        ros(nei)= rho

        p=   PA(1,j)
        u=   PA(2,j)
        v=   PA(3,j)
        tem= PA(4,j)
        y(1:Nvar-Nflow)=PA(nflow+1:Nvar,j); y(Nvar-Nflow+1)= 1. - sum(y(1:Nvar-Nflow))
        call ComputeGasParameter(y, tem, Cp, gama, Rcpcv)
        rho= p/(Rcpcv*tem)
        romax= max(romax,rho)
      enddo
      
      sum1=0.; 
      do nei=1, Neigh(i)   !! 3
        sum1=sum1+ abs(ros(nei)-ros(0))
      enddo
      Smooth(i)=sum1/(Neigh(i)*romax*vol(i)*rc_detector)
    enddo
    
    return
end subroutine
    