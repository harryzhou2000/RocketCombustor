!C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ReadMesh   ! MESH FLOW
    use mainvar
    use Parallel_Var
    use ConnectFace
    implicit none
    integer i,k,ic,n1,n2,np,ib,icon,ibound_max
    type(Connect), pointer :: tCon
    
    filename= trim("grid"//trim(strProc)//".in")
    OPEN(5,FILE=filename,STATUS='OLD')
!! Read mesh on cell-based format
    read(5,*) NE,NE_Vir, NNT,NNT_Vir,NF,NF_Vir
    write(*,*) NE, " ", NE_Vir, "NE" 
    
    call allocate_cellandnode_arrary
    
    call allocate_face_arrary
    
    do i=1,NE_Vir
      read(5,*) kface(i), (N(k,i),k=1,kface(i))
    enddo
    
    call allocate_nodetopology_arrary
      
    call allocate_stencil_arrary
    
!! read coordinates
    READ(5,*)(COOR(1,NP),COOR(2,NP),IsCommonBoundaryNode(1,NP),IsCommonBoundaryNode(2,NP),NP=1,NNT_Vir)
    Coor= Coor /unit
    
    !! reset N(4,i)
    do i=1,NE_Vir
    if(kface(i)==3) N(4,i)= N(3,i)
    enddo
    
    !! read connectivity info
    read(5,*) NConE
    allocate(ConE(NConE))
    do icon=1, NConE
      tCon=> ConE(icon)
      read(5,*) tCon%tarPar,tCon%tarCon,tCon%NListR,tCon%NListS
      allocate(tCon%ListR(2,tCon%NListR))
      allocate(tCon%ListS(2,tCon%NListS))
      read(5,*) ( tCon%ListR(1,ic),tCon%ListR(2,ic),ic=1,tCon%NListR )
      read(5,*) ( tCon%ListS(1,ic),tCon%ListS(2,ic),ic=1,tCon%NListS )
      write(*,*) NConE, "Here" 
    enddo
    
!! read boundary infomation
    read(5,*) NBoundary,ibound_max
    
    if(NBoundary>0) call allocate_boundary_arrary(ibound_max)
    
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

    !! calculate the centroid and volume
    call VolumeAndCenter

    call find_monitor_point
    
    !! output grids
    filename= trim("grid"//trim(strProc)//".plt")
    open(1,FILE=filename)
    WRITE(1,*)'VARIABLES = "X", "Y"'
    WRITE(1,*)'ZONE N=',NNT,',E=',NE,', F=FEPOINT, ET=QUADRILATERAL'
    DO I=1,NNT
    WRITE(1,*) COOR(1,I),COOR(2,I)
    ENDDO

    DO I=1,NE
    WRITE(1,*) (N(k,i),k=1,4)
    ENDDO
    CLOSE(1)
      
    return
end subroutine



!C   *****************************************************************
      subroutine change
      use mainvar
      real*8 :: xv1(2),xv2(2),xv3(2),xv4(2)
      integer:: numofface,numofface_vir
      
      !! create faces
      Neighbor(:,:)= -100

      i=0  !! i �����Ǳߵ�������
      do j=1,NE_Vir
      if( mod(j,10000)==0 ) write(*,*) 'cell j:',j
        do kf=1,KFace(j)
         n1= N(kf,j)
         n2= N(mod(kf,KFace(j))+1,j)
         call face(n1,n2,i,j )  !�Ե�kf����
        enddo
        if( j==NE ) numofface= i
      enddo
      numofface_vir=i
      
      if(numofface /= NF) then
        write(*,*) 'error in number of face'
        stop
      endif
      
      if(numofface_vir /= NF_vir) then
        write(*,*) 'error in number of virtual face'
        stop
      endif
      
      write(*,*) 'number of faces:', NF,NF_Vir
      
      return
end subroutine

!c--------------------------------------------------------------------------
subroutine face(n1,n2,i,j)
      use mainvar
      
      IndexFace= line(n2,n1)

      if( IndexFace .eq. 0) then         
        i=i+1
        neighbor(1,i)=j     !! ��Ϊj��Ԫ����ʱ�밲�Ŷ��㣬then cell j is side right cell
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



!!---------------------------------------------
!! volume and gravity of cell
!!--------------------------------------------
subroutine VolumeAndCenter
    use mainvar
    implicit none
    include 'mpif.h'
    double precision:: xv1(2),xv2(2),xv3(2),xv4(2),dxdy(2)
    double precision:: sumvol,sumvol_sum
    integer :: i,inode,jnode,idirect,ierr
    double precision,external:: GaussInteg_Quad
    
          
    do i=1,NE_vir
      
      xv1(:)= Coor(:,n(1,i))
      xv2(:)= Coor(:,n(2,i))
      xv3(:)= Coor(:,n(3,i))
      xv4(:)= Coor(:,n(4,i))
       

      vol(i)= GaussInteg_Quad(fun, xv1,xv2,xv3,xv4)
      ! write(*,*) xv1, xv2, xv3, xv4, " Vol Coords: ", vol(i)
      cellxy(:,i)= (/ GaussInteg_Quad(funx, xv1,xv2,xv3,xv4)/vol(i), &
                     GaussInteg_Quad(funy, xv1,xv2,xv3,xv4)/vol(i) /)

       
      if( vol(i)<0. )then
        write(*,*) 'negative volume', i, vol(i)
        stop
      endif
      
!!======================================================
!!  reference length scale for non-dimension
!!------------------------------------------------------
      
       reflen(:,i)=0.

       do inode=1,Kface(i)-1
       do jnode=inode+1,Kface(i)
       do idirect=1,2
       reflen(idirect,i)=max(0.5*dabs(Coor(idirect,n(inode,i))-Coor(idirect,n(jnode,i))),reflen(idirect,i))
       enddo
       enddo
       enddo
       
       !do inode=1,Kface(i)
       !do idirect=1,2
       !reflen(idirect,i)=max(abs(Coor(idirect,n(inode,i))-cellxy(idirect,i)),reflen(idirect,i))
       !enddo
       !enddo
      
    enddo

    call Connect_CellCenter  !! check connected cells

    sumvol=0.
    do i=1,NE
      sumvol= sumvol + vol(i)
    enddo

    call MPI_AllReduce(sumvol,sumvol_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

    total_volume= sumvol_sum

    write(*,*) 'Total volume of the mesh is ', total_volume
    
    
contains
    real function fun(xvec)
    real :: xvec(2)
    fun= 1.
    end function
    
    real function funx(xvec)
    real :: xvec(2)
    funx= xvec(1)
    end function
    real function funy(xvec)
    real :: xvec(2)
    funy= xvec(2)
    end function
    
end subroutine
