subroutine compute_reconstruction_matrix
    use mainvar
    implicit none
    integer :: NTemp(NE_Vir)
    integer :: i,j,k,nei,kf,ic,iface,iorder,NL,NR,Ninter,Nbound,Log
    real*8  :: Offset(2),xv1(2),xv2(2),xv3(2),xv4(2)
    real*8,allocatable :: A1(:,:), A1inv(:,:) 
    real*8  :: distance 
    real*8,external :: GaussInteg_Quad
    
!! NeighF(3,NE)
    Ntemp(:)=0
    do i=1, NF_Vir
      NL = Neighbor(1,i)
      NR = Neighbor(2,i)
      
      if(NL<=0)stop
      Ntemp(NL)=Ntemp(NL)+1
      NeighF(Ntemp(NL),NL)= i
      if(NR>0 .and. Fproperty(i)/=100)then
        Ntemp(NR)=Ntemp(NR)+1
        NeighF(Ntemp(NR),NR)= i
      endif
    enddo
    
    
!! first neighbor cells: faces neighbored
    neigh(:) = 0
    Ninter=0
    NBound=0
    do i=1, NE
      do k=1,KFace(i)  !! 3
        iface = NeighF(k,i)    !! Î´ÓÐ¶¨Òå
        do kf=1,2
        ic= Neighbor(kf,iface)
        if( ic/=i .and. ic>0 )then
          Neigh(i)= Neigh(i)+1
          NeighC(Neigh(i),i)= ic
          
          Ofs(1:2,Neigh(i),i)= SideOfs(1:2,iface)
        endif
        enddo
      enddo
      
      if(neigh(i)==kface(i)) then
       Ninter=Ninter+1
      else
       Nbound=Nbound+1
      endif
    enddo

    do i=1, NE
      NeighC (0,i)= i
      Ofs(1:2,0,i)= 0.
    enddo
    
    if(Nbound+Ninter/=NE) then
    write(*,*) 'error in boundcell number'
    stop
    endif 

!!------------------------------------------
    write(*,*) 'reconstruction stencil,min,mincell/max,maxcell:',minval(Neigh(1:NE)),minloc(Neigh(1:NE)),&
                maxval(Neigh(1:NE)),maxloc(Neigh(1:NE))


    !! construct ghost cells
    call construct_ghostcells
    
    
!! moments in basis
    do i=1, NE_Vir
      
      xv1=Coor(:,N(1,i)) 
      xv2=Coor(:,N(2,i)) 
      xv3=Coor(:,N(3,i)) 
      xv4=Coor(:,N(4,i))

      do iorder=1,MaxOrder
      XYH(iorder,i)= GaussInteg_Quad(fun, xv1,xv2,xv3,xv4)/vol(i)
      enddo

      if( abs(XYH(1,i))+abs(XYH(2,i))>1.e-8 )then
        write(*,*) 'error in XYH,', i, abs(XYH(1,i)),abs(XYH(2,i)),vol(i)
        stop
      endif
    enddo
    
!!======================================================
!!  Form matrices for least-squares reconstruction
!!------------------------------------------------------
    
    !! construct least-square matrix
    write(*,*) 'constructing least-squares reconstruction matrix ...'
    
    do i=1, NE

      allocate(A1   (kface(i),MaxOrder)); A1=0.
      allocate(A1inv(MaxOrder,kface(i))); A1inv=0.
       
      do k=1, Neigh(i)
        Offset(1:2)= Ofs(1:2,k,i) !! offset

        j= NeighC(k,i)

        xv1=Coor(:,N(1,j))-Offset(:) 
        xv2=Coor(:,N(2,j))-Offset(:) 
        xv3=Coor(:,N(3,j))-Offset(:) 
        xv4=Coor(:,N(4,j))-Offset(:) 

        distance= sqrt( (CellXY(1,j)-Offset(1) -CellXY(1,i))**2+ &
                      (CellXY(2,j)-Offset(2) -CellXY(2,i))**2 )
        wdis(k,i)= distance**(-1)

        do iorder=1, MaxOrder
          A1(k,iorder)=wdis(k,i)*( GaussInteg_Quad( fun, xv1,xv2,xv3,xv4)/vol(j) - XYH(iorder,i) )
        enddo
      enddo

      if(Num_ghostcell(i)>0) then
        do k=1, Num_ghostcell(i)
          j= Ghostcell(k,i)  

          xv1=Coor_Ghostnode(:,1,j)
          xv2=Coor_Ghostnode(:,2,j)
          xv3=Coor_Ghostnode(:,3,j)
          xv4=Coor_Ghostnode(:,4,j)

          distance= sqrt( (CellXY_Ghostcell(1,j)-CellXY(1,i))**2+ &
                          (CellXY_Ghostcell(2,j)-CellXY(2,i))**2 )
          wdis(k+neigh(i),i)= distance**(-1)

          do iorder=1, MaxOrder  
            A1(k+neigh(i),iorder)=wdis(k+neigh(i),i)*( GaussInteg_Quad( fun, xv1,xv2,xv3,xv4)/vol_Ghostcell(j) - XYH(iorder,i) )
          enddo
        enddo
      endif
      
      ! call BGINV(kface(i),MaxOrder,A1,A1inv, Log, 1.e-15)
      ! if( Log==1 ) write(*,*) " error in bginv!! "

      call pseudo_inverse_lapack(A1,kface(i),MaxOrder,A1inv)
      
      AMLS(1:MaxOrder,1:kface(i),i) = A1Inv(1:MaxOrder,1:kface(i))
      
      
      deallocate(A1)
      deallocate(A1inv)
      
    enddo
    
contains
    real*8 function fun(xy)
    implicit none
    real*8 :: xy(2),x,y
      x=(xy(1)-CellXY(1,i))/reflen(1,i)
      y=(xy(2)-CellXY(2,i))/reflen(2,i)
    select case(iorder)
    case(1)
    fun= x
    case(2)
    fun= y
    case default
      write(*,*) 'error in iorder!',iorder,i; stop
    end select
    return
    end function
    
end subroutine


logical function contain(in,iarr,n)
implicit none
integer :: in,n,iarr(n),i
do i=1,n
if( in==iarr(i) )then
  contain= .true.
  return
endif
enddo
contain=.false.
end function
    

