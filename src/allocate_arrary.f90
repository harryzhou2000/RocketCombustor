!!========================================
!! Allocate flow solver related arraries
!!----------------------------------------
subroutine allocate_cellandnode_arrary
    use mainvar
    implicit none

 !! mesh size info
    integer  NCell,NEU
    
    Ncell= NE_Vir
    NEU  = NNT_Vir
    
    !! Mesh topological relation

    !! cell
    allocate( N(4,Ncell) )
    allocate( NeighF(4,Ncell) ) 
    allocate( KFace(Ncell) )
    
    !! node coordinates
    allocate( COOR(2,NEU) ) ; coor=0.
    !! cell 
    allocate( VOL(Ncell) ) ; vol=0.
    allocate( CELLXY(2,Ncell) ); cellxy=0.
    allocate( reflen(2,Ncell) );reflen=0.
    
    !! physical data, reconstruction coefficients,flux
    allocate( WA(Nvar,Ncell) ) ; WA=0.
    allocate( PA(Nvar,Ncell) ) ; PA=0.
    allocate( WP(Nvar,Ncell) ) ; WP=0.
    allocate( WP2(Nvar,Ncell) ) ; WP2=0.
    allocate( GradC(MaxOrder,NVar,NCell) ) ; gradC=0.
    allocate( grad_tmp(MaxOrder,NVar,NCell) ) ;grad_tmp=0.
    allocate( hWA(Nvar,Ncell) ) ; hWA=0.
    allocate( hWAP(InLoop,Nvar,Ncell) ) ; hWAP=0.
    allocate( XYH(MaxOrder,NCell) ) ;xyh=0.
    allocate( Smooth(Ncell) );smooth=0.
    allocate( Deltaeig(Ncell) ); Deltaeig=0.

    !! combustion reaction rate
    allocate(production_rate(NVar-NFlow,Ncell)); production_rate=0. 
    allocate(production_rate_sum(Ncell)); production_rate_sum=0. 
    allocate(flux_difference(Ncell)); flux_difference=0.
    allocate(flux_difference_nosource(Ncell)); flux_difference_nosource=0.
    allocate(sum_of_mass_fraction(Ncell)); sum_of_mass_fraction=0.
    allocate(phi_cell(Nvar,Ncell)); phi_cell=0.
   
    !! local time step
    allocate( step(Ncell) ) ; step=0.
    
    !! for implicit RK 
    allocate( fixed(Nvar,Ncell) ) ; fixed=0.

    !! commen boundary information
    allocate(IsCommonBoundaryNode(2,NEU)); IsCommonBoundaryNode=0

    return
end subroutine
    
subroutine allocate_stencil_arrary
    use mainvar
    implicit none
    integer  Ncell
    
    Ncell=NE
    
    MaxNei=maxval(Kface(:))
    
    !! Mesh topological relation
    !! cell
    allocate( NeighC (0:MaxNei,Ncell) ); NeighC=0
    allocate( Neigh (Ncell) ) ; Neigh=0

    allocate( wdis(MaxNei,Ncell) ); wdis=1.

    allocate( Ofs(2,0:MaxNei,NCell) ); Ofs=0.
    
    allocate( AMLS(MaxOrder,MaxNei,1:NE)); AMLS=0.

return
end subroutine
    
subroutine allocate_nodetopology_arrary
    use mainvar
    implicit none
    integer i,k,nn,max_nodeface
    
      allocate( NNodeNC(NNT_vir) )
    
      NNodeNC(:) =0
      do i=1, NE_Vir
      do k=1, KFace(i)
      nn= N(k,i)
      NNodeNC(nn)= NNodeNC(nn)+1
      enddo
      enddo
      
      allocate( NodeNC(maxval(NNodeNC(:)),NNT_Vir) )
      
      NNodeNC(:) =0
      do i=1, NE_Vir
      do k=1, KFace(i)
      nn= N(k,i)
      NNodeNC(nn)= NNodeNC(nn)+1
      NodeNC( NNodeNC(nn),nn ) = i
      enddo
      enddo
      
      max_nodeface=2*maxval(NNodeNC(:))
      
      allocate( NodeFace(max_nodeface,NNT_vir))
      allocate( NumNFace(NNT_vir))             ; NumNFace=0

return
end subroutine

subroutine allocate_face_arrary
    use mainvar
    implicit none
    
    !! Mesh topological relation
   
    !!face
    allocate( FNode(2,NF_vir) )
    allocate( NEIGHBOR(2,NF_vir) )
    allocate( FProperty(NF_vir) ); FProperty=0

    !! face
    allocate( VECX(2,NF_vir) ) ; VECX=0.
    allocate( SideOfs(2,NF_vir) );SideOfs=0.
   

    !! spectral radius of cell interfaces
    allocate( alpaF(NF_vir) ) ; alpaF=0.

    !! common boundary face
    allocate(IsCommonBoundaryFace(NF)); IsCommonBoundaryFace=0

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
    
    !! periodic boundaries
    allocate(Con(2,NBoundary))
    allocate( dxyz(2,NBoundary) );dxyz=0.

return
end subroutine