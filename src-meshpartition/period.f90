!!======================================================
!!      periodic cell and boundary at every time step 
!!=====================================================

module ConnectFace
    !! connectivity faces
    !! all boundaries
    integer,save :: NBoundary
    integer,allocatable :: iBound(:),iType(:),BoundN(:,:,:),BoundF(:,:),BoundFT(:,:)
    integer,save :: NCon
    integer,allocatable:: Con(:,:)
    real*8,allocatable:: dxyz(:,:)
end module    

subroutine CreatePeriodConnect
    use mainvar
    use ConnectFace
    integer,allocatable ::  NNodeNF(:),NodeNF(:,:)
    real :: ddd(2),xfs(2),xft(2)
    logical,external :: contain
    real,external :: sumup
    
    allocate( NNodeNF(NNT) )
    
    NNodeNF(:) =0
    do i=1, NF
    do k=1, 2
      nn= FNode(k,i)
      NNodeNF(nn)= NNodeNF(nn)+1
    enddo
    enddo
    
    allocate( NodeNF(maxval(NNodeNF(:)),NNT) )
    
    NNodeNF(:) =0
    do i=1, NF
    do k=1, 2
      nn= FNode(k,i)
      NNodeNF(nn)= NNodeNF(nn)+1
      NodeNF(NNodeNF(nn),nn)=i
    enddo
    enddo
    
    
!! find boundary face
    do ib=1,NBoundary
    do np=1,iBound(ib)
      n1=BoundN(1,np,ib)
      n2=BoundN(2,np,ib)
      !! 找边界面
      do ic1=1,NNodeNF(n1)
        iface= NodeNF(ic1,n1)
        if( Contain(n2, FNode(:,iface),2) )then
        BoundF(np,ib)= iface
        exit
        endif
      enddo
      if( ic1>NNodeNF(n1) )then
        write(*,*) 'not find boundary face',ib,np;
        stop
      endif
    enddo
    enddo
    
    if( NCon==0 ) return
    
!! create face correlation,特别是最后的Neighbor(2,iface),为后面的找模板提供方便
    do icon=1,NCon
      is= Con(1,icon)
      it= Con(2,icon)
      do ib=1,iBound(is)
        
        !! 先将面对应的目标单元找到
        n1= BoundN(1,ib,is)
        n2= BoundN(2,ib,is)
        xfs(:)= (Coor(:,n1)+ Coor(:,n2))*0.5
        
        do itb=1,iBound(it)
          n1= BoundN(1,itb,it)
          n2= BoundN(2,itb,it)
          sav= sqrt((Coor(1,n1)-Coor(1,n2))**2 + &
                    (Coor(2,n1)-Coor(2,n2))**2)
          xft(:)= ( Coor(:,n1)+Coor(:,n2) )*0.5
          !! n1 ivs1,  n2 ivs2 how to judge they are periodic
          ddd(:)= xft(:)- xfs(:)
          ! if( abs(ddd(1)-dxyz(1,icon))+ abs(ddd(2)-dxyz(2,icon))< 1.0e-7 )  exit
          if(sqrt((ddd(1)-dxyz(1,icon))**2 + &
                  (ddd(2)-dxyz(2,icon))**2) < 0.1*sav) exit
        enddo
        if(itb>iBound(it))then
          write(*,*) 'not find target.'
          stop
        endif
        BoundFT(ib,is)= BoundF(itb,it)  !! create connect face relation
        iFace= BoundF(ib, is)
        itFace=BoundF(itb,it)
        Neighbor(2,iface)= Neighbor(1,itFace)
        
        FProperty(iface)= 100
        SideOfs(:,iface)= dxyz(:,icon) !!
      enddo
    enddo
    
    deallocate( NNodeNF)
    deallocate( NodeNF )
    
end subroutine
