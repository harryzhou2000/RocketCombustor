!C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INPUT   ! MESH FLOW
    use mainvar
    use ConnectFace
    implicit none
    integer i,k,n1,n2,np,ib,icon,ibound_max
    
    OPEN(5,FILE='grid.in',STATUS='old')
! Read mesh on cell-based format
    read(5,*) NE,NNT
    
    call allocate_cellandnode_arrary
    
    do i=1,NE
      read(5,*) kface(i), (N(k,i),k=1,kface(i))
    enddo
    
    call allocate_nodetopology_arrary
      
    call allocate_stencil_arrary
    
!! read coordinates
    ! READ(5,*)(COOR(1,NP),COOR(2,NP),NP=1,NNT)
    ! READ(5,*)(COOR(1,NP),COOR(2,NP),IsCommenBoundaryNode(1,NP),IsCommenBoundaryNode(2,NP), NP=1,NNT)
    do i=1,NNT
      read(5,*) COOR(1,i),COOR(2,i),IsCommenBoundaryNode(1,i),IsCommenBoundaryNode(2,i)
    enddo
    Coor= Coor /unit
    
    !! reset N(4,i)
    do i=1,NE
    if(kface(i)==3) N(4,i)= N(3,i)
    enddo
    
!! read boundary infomation
    read(5,*) NBoundary,ibound_max
    
    call allocate_boundary_arrary(ibound_max)
    
    do ib=1, NBoundary
      read(5,*) itype(ib),ibound(ib)
      read(5,*) (BoundN(1,np,ib),BoundN(2,np,ib),np=1,ibound(ib))
    enddo
    
    read(5,*) NCon
    do icon=1,NCon
      read(5,*) (Con(k,icon),k=1,2), (dxyz(k,icon),k=1,2)
      iType( Con(1,icon) )= 100
      iType( Con(2,icon) )= 100
    enddo
    
    CLOSE(5)
    
    call allocate_face_arrary

    !! create faces
    call change
    
    !! init periodic connectivity
    call CreatePeriodConnect
    
  !! set face property,  inner/boundary type
	do i=1, NF
      n1= FNode(1,i); n2= FNode(2,i)
	  if( Neighbor(2,i)<0 )then !! boundary
	    do ib=1, NBoundary
        do np=1, ibound(ib)
          if( (n1==BoundN(1,np,ib) .and. n2==BoundN(2,np,ib))  .or. &
              (n2==BoundN(1,np,ib) .and. n1==BoundN(2,np,ib)) )then
            Neighbor(2,i)= -iType(ib)
            exit
          endif
        enddo
          if(np<=ibound(ib))exit
        enddo
        if( ib>NBoundary )then
          write(*,*) 'not find boundary of this type'; stop
        endif
	  endif
	enddo
    
	!! output grids
	  OPEN(1,FILE='grid.PLT') 
      WRITE(1,*)'VARIABLES = "X", "Y"'
      WRITE(1,*)'ZONE N=',NNT,',E=',NE,', F=FEPOINT, ET=QUADRILATERAL'
      DO I=1,NNT
      WRITE(1,*) COOR(1,I),COOR(2,I)
      ENDDO

      DO I=1,NE
      WRITE(1,*) (N(k,i),k=1,4)
      ENDDO
      CLOSE(1)
      
end subroutine



!C   *****************************************************************
      subroutine change
      use mainvar
      real*8 :: xv1(2),xv2(2),xv3(2),xv4(2)
      integer:: numofface,numofface_vir
      
      !! create faces
      Neighbor(:,:)= -100

      i=0  !! i 用来记边的累增量
      do j=1,NE
      if( mod(j,10000)==0 ) write(*,*) 'cell j:',j
        do kf=1,KFace(j)
         n1= N(kf,j)
         n2= N(mod(kf,KFace(j))+1,j)
         call face(n1,n2,i,j )  !对第kf条边
        enddo
      enddo
      
      numofface= i
      
      if(numofface /= NF) then
        write(*,*) 'error in number of face'
        stop
      endif
      
      write(*,*) 'nf:', NF
      
      return
end subroutine

!c--------------------------------------------------------------------------
subroutine face(n1,n2,i,j)
      use mainvar
      
      IndexFace= line(n2,n1)

      if( IndexFace .eq. 0) then         
        i=i+1
        neighbor(1,i)=j     !! 因为j单元是逆时针安排顶点，then cell j is side right cell
        vecx(1,i)=  coor(2,n2)-coor(2,n1)  !! vector of side i
        vecx(2,i)=-(coor(1,n2)-coor(1,n1))

	   FNode(1,i)= n1
	   FNode(2,i)= n2
       
       !! add face index to node
       NumNFace(n1)= NumNFace(n1)+1
       NumNFace(n2)= NumNFace(n2)+1
       if( NumNFace(n1)>50 .or. NumNFace(n2)>50 )then
         write(*,*) 'error in NumNFace'
         stop
       endif
       NodeFace(NumNFace(n1),n1)= i
       NodeFace(NumNFace(n2),n2)= i
       
      else
        neighbor(2,IndexFace)=j
      endif
    
    
	!! line(n2,n1), line(n1,n2) is inverse direction vector
contains
   integer function Line(n1,n2)
   do inum=1,NumNFace(n1)
     iface= NodeFace(inum,n1)
     if( (FNode(1,iface)==n1 .and. FNode(2,iface)==n2) .or. &
         (FNode(1,iface)==n2 .and. FNode(2,iface)==n1) )then
       Line= iface
       return
     endif
   enddo
   Line=0
   return
   end function

end subroutine
