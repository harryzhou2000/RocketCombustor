subroutine allocate_cellandnode_arrary
    use mainvar
    implicit none

 !! mesh size info
    integer  NCell,NEU
    
    Ncell= NE
    NEU  = NNT
    
    !! Mesh topological relation

    !! cell
    allocate( N(4,Ncell) )
    allocate( NeighF(4,Ncell) ) 
    allocate( KFace(Ncell) )
    
    !! node coordinates
    allocate( COOR(2,NEU) )
    !! cell 
    allocate( VOL(Ncell) )
    allocate( CELLXY(2,Ncell) )
    
    !! curve and wall cell index
    allocate( IsCurveCell(Ncell) ); IsCurvecell=-1
    allocate( IfWallcell(Ncell) ) ; IfWallcell=-1
    
    !! bound info
    allocate( IsBoundCell(Ncell) ); IsBoundCell(:)=-1
    allocate( itib(Ncell) )
    allocate( inter(NCell))
    allocate( bound(NCell))


    !! commen boundary information
    allocate(IsCommenBoundaryNode(2,NEU));IsCommenBoundaryNode=0

    return
end subroutine
    
subroutine allocate_stencil_arrary
    use mainvar
    implicit none
    integer  Ncell
    
    Ncell=NE
    
    MaxNei=maxval(Kface(:))
    MaxNei_bound=MaxNei-1
    
    !! Mesh topological relation
    !! cell
    allocate( NeighC (0:MaxNei,Ncell) )
    allocate( Neigh (Ncell) )

    allocate( Ofs(2,0:MaxNei,NCell) );Ofs=0.

return
end subroutine
    
subroutine allocate_nodetopology_arrary
    use mainvar
    implicit none
    integer i,k,nn,max_nodeface
    
      allocate( NNodeNC(NNT) )
    
      NNodeNC(:) =0
      do i=1, NE
      do k=1, KFace(i)
      nn= N(k,i)
      NNodeNC(nn)= NNodeNC(nn)+1
      enddo
      enddo
      
      allocate( NodeNC(maxval(NNodeNC(:)),NNT) )
      
      NNodeNC(:) =0
      do i=1, NE
      do k=1, KFace(i)
      nn= N(k,i)
      NNodeNC(nn)= NNodeNC(nn)+1
      NodeNC( NNodeNC(nn),nn ) = i
      enddo
      enddo
      
      max_nodeface=2*maxval(NNodeNC(:))
      
      allocate( NodeFace(max_nodeface,NNT))
      allocate( NumNFace(NNT))             ; NumNFace=0

return
end subroutine

subroutine allocate_face_arrary
    use mainvar
    use ConnectFace
    implicit none

 !! mesh size info
    integer  Ntotalface,NumofFace,NumofBoundFace,Nface
    integer  i,ib,ifacecount
    
    Ntotalface=0
    do i=1,NE
      select case( Kface(i) )
      case(3)  !! triangle
        Ntotalface=Ntotalface+3
      case(4)  !! quad
        Ntotalface=Ntotalface+4
      case default
        write(*,*) 'error. no such type of cell'
        stop
      end select
    enddo
    
    NumofBoundFace=0
    do ib=1, NBoundary
	  NumofBoundFace=NumofBoundFace+ibound(ib)
    enddo
    
    NumofFace=(Ntotalface+NumofBoundFace)/2
    NF=NumofFace
    
    !! Mesh topological relation
   
    !!face
    allocate( FNode(2,NF) )
    allocate( NEIGHBOR(2,NF) )
    allocate( FProperty(NF) )

    !! face
    allocate( VECX(2,NF) )
    allocate( SideOfs(2,NF) );SideOfs=0.
   
    !! curve face index
    allocate( IsCurveFace(NF) ); IsCurveFace=-1

return
end subroutine
    
subroutine allocate_boundary_arrary(ibound_max)
    use ConnectFace
    implicit none
    integer,intent(in)::  ibound_max

    allocate(iBound(NBoundary))
    allocate(iType(NBoundary))
    allocate(BoundN(2,ibound_max,NBoundary))
    allocate( BoundF(ibound_max,NBoundary))
    allocate(BoundFT(ibound_max,NBoundary))
    
    !! 注意，连接面数组，给多了
    allocate(Con(2,NBoundary))
    allocate( dxyz(2,NBoundary) );dxyz=0.

return
end subroutine