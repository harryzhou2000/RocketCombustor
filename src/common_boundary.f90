!!======================================
!! output the common boundaries data
!!--------------------------------------
subroutine save_common_boundary
  use mainvar
  use Parallel_Var
  implicit none
  integer,parameter:: IOfile1= 115, IOfile2=151
  double precision:: xy(2),dxyL(2),dxyR(2),vvh(MaxOrder),sol(Nvar)
  integer:: i,j,k,n1,n2,NL,NR,Nswap
  integer,save :: isave=0
  character*50 :: filename1, filename2

  if(IfRightDomain) return
  if(IfSteady) return
  if(ttime < save_start_time - 0.1*dt .or. ttime > save_end_time + 0.1*dt) return

  
  !! output u_l
  filename1= trim("L_solution_"//trim(strProc)//".dat")
  filename2= trim("R_solution_"//trim(strProc)//".dat")
  if(isave>0) then
    open(IOfile1, file=filename1, status="old", position="append", action="write")
    open(IOfile2, file=filename2, status="old", position="append", action="write")
  else
    open(IOfile1, file=filename1)
    open(IOfile2, file=filename2)
  endif

  ! call reconstruction

  do i=1,NF
    if(IsCommonBoundaryFace(i)>0) then
      n1= FNode(1,i)
      n2= FNode(2,i)
      xy(:)= 0.5*(Coor(:,n1)+Coor(:,n2))
      NL= Neighbor(1,i)
      NR= Neighbor(2,i)

      if(CellXY(1,NL)>CellXY(1,NR)) then
        Nswap= NL
        NL= NR
        NR= Nswap
      endif

      if(NL>NE) cycle

      !! Left and Right state
      dxyL(:)= xy(:)- CellXY(:,NL)
      call fbaseAll(dxyL,vvh,NL)

      do k=1,Nvar
        sol(k)= PA(k,NL) + sum(gradC(:,k,NL)*vvh(:))
      enddo

      write(IOfile1,*) ttime,sol(:)

      dxyR(:)= xy(:)- CellXY(:,NR)
      call fbaseAll(dxyR,vvh,NR)

      do k=1,Nvar
        sol(k)= PA(k,NR) + sum(gradC(:,k,NR)*vvh(:))
      enddo

      write(IOfile2,*) ttime,sol(:)

    endif
  enddo

  close(IOfile1)
  close(IOfile2)

  isave=isave + 1

  return
end subroutine


!!===========================================
!! output the common boundaries information
!!-------------------------------------------
subroutine save_common_boundary_info
  use mainvar
  use Parallel_Var
  implicit none
  integer,parameter:: IOfile= 115, num_parameters=0
  double precision:: xy(2),dxyL(2),rof,uf,vf,pf,tf,vvh(MaxOrder)
  integer:: i,j,k,n1,n2,NL,NR,num_boundary_face

  if(IfRightDomain) return

  open(unit=IOfile,file="info_solution_"//trim(strProc)//".dat")

  write(IOfile,*) num_parameters

  ! write(IOfile,*) num_common_boundary

  num_boundary_face=0
  do i=1,NF
    if(IsCommonBoundaryFace(i)>0) then
      n1= FNode(1,i)
      n2= FNode(2,i)
      xy(:)= 0.5*(Coor(:,n1)+Coor(:,n2))
      NL= Neighbor(1,i)
      NR= Neighbor(2,i)

      if(CellXY(1,NL)>CellXY(1,NR)) NL= NR

      if(NL>NE) cycle

      num_boundary_face=num_boundary_face+1

    endif
  enddo

  write(IOfile,*) num_boundary_face

  write(*,*) 'there are ', num_boundary_face, 'common boundary elements in processor ', myProc

  do i=1,NF
    if(IsCommonBoundaryFace(i)>0) then
      n1= FNode(1,i)
      n2= FNode(2,i)
      xy(:)= 0.5*(Coor(:,n1)+Coor(:,n2))
      NL= Neighbor(1,i)
      NR= Neighbor(2,i)

      if(CellXY(1,NL)>CellXY(1,NR)) NL= NR

      if(NL>NE) cycle

      write(IOfile,*) IsCommonBoundaryFace(i),',',xy(1),',',xy(2)

    endif
  enddo

  close(IOfile)

  return
end subroutine


!! detect the common boundaries
subroutine detect_common_boundary
  use mainvar
  use Parallel_Var
  implicit none
  integer, parameter :: max_num_common_boundary=100
  integer:: common_boundary(max_num_common_boundary)
  integer:: i,j,k,n1,n2

  num_common_boundary=0
  common_boundary=0

  do i=1,NF
    n1= FNode(1,i)
    n2= FNode(2,i)

    if(IsCommonBoundaryNode(1,n1)*IsCommonBoundaryNode(2,n1)>0) then  !! n1 belongs to two common boundaries
      if(IsCommonBoundaryNode(1,n2)>0 .and. IsCommonBoundaryNode(2,n2)==0) then
        IsCommonBoundaryFace(i)= IsCommonBoundaryNode(1,n2)
      else
        write(*,*) 'error in the setting of common boundary point:',n2,'in processor:', MyProc
        call MPI_STOP
      endif

    elseif(IsCommonBoundaryNode(1,n1)+IsCommonBoundaryNode(2,n1)>0) then !! n1 just belongs to one common boundaries
      if(IsCommonBoundaryNode(1,n2)+IsCommonBoundaryNode(2,n2)>0) then   !! n2 also belongs to common boundary
        if(any(IsCommonBoundaryNode(1,n1)==IsCommonBoundaryNode(:,n2))) then 
          IsCommonBoundaryFace(i)= IsCommonBoundaryNode(1,n1)
        else
          write(*,*) 'error in the setting of common boundary point:',n2,'in processor:', MyProc
          call MPI_STOP
        endif
      endif

    endif

    if(IsCommonBoundaryFace(i)>0) then
      if(.not. any(IsCommonBoundaryFace(i)==common_boundary)) then
        num_common_boundary = num_common_boundary +1
        common_boundary(num_common_boundary)= IsCommonBoundaryFace(i)
      endif
    endif

  enddo


  write(*,*) 'there are', num_common_boundary, 'common boundaries in processor', MyProc

  return
end subroutine


subroutine arrange_common_boundary_cell_index
  use mainvar
  use Parallel_Var
  implicit none
  include 'mpif.h'
  integer, parameter :: IOfile= 114
  double precision,parameter :: epsilon=1.0e-10
  double precision:: xy(2)
  integer:: ipar,ic,id,num_common_boundary_cells,icount_common_boundary_cell,icount,&
            NL,NR,n1,n2,icount_sum,ierr,iface

  integer,allocatable:: find_common_boundary_cell(:),id_face(:,:)
  double precision,allocatable:: coor_common_boundary_cell(:,:)
  character*50  str

  allocate(find_common_boundary_cell(NF));find_common_boundary_cell=0
  allocate(index_common_boundary(2,NF));index_common_boundary=0
  ! allocate(common_boundary_xy(2,NF));common_boundary_xy=0.

  icount_common_boundary_cell= 0
  icount= 0
  do ipar=1,NPreviousPar
    if(ipar<10)then
      write(str,"(i1)") ipar
    elseif(ipar>=10 .and. ipar<100) then
      write(str,"(i2)") ipar
    elseif(ipar>=100 .and. ipar<1000) then
      write(str,"(i3)") ipar
    elseif(ipar>=1000 .and. ipar<10000) then
      write(str,"(i4)") ipar
    else
      write(*,*) 'Error: processor index', ipar, 'larger than 10000 when looking for the common boundary information file'
      call MPI_STOP
    endif

    open(unit=IOfile,file="face_index_"//trim(str)//".dat",status="old",action="read")
    read(IOfile,*) num_common_boundary_cells

    if(num_common_boundary_cells==0) then
      close(IOfile)
    else
    
      allocate(id_face(2,num_common_boundary_cells));id_face=0
      allocate(coor_common_boundary_cell(2,num_common_boundary_cells));coor_common_boundary_cell=0.
      do ic=1,num_common_boundary_cells
        read(IOfile,*) id_face(1,ic),id_face(2,ic), coor_common_boundary_cell(1,ic),coor_common_boundary_cell(2,ic)
      enddo
      close(IOfile)

      do iface=1,NF
        NL= neighbor(1,iface)
        NR= neighbor(2,iface)
        n1= FNode(1,iface)
        n2= FNode(2,iface)
        xy(:)= 0.5*(Coor(:,n1)+Coor(:,n2))
        if(-NR>1000) then
          do ic=1,num_common_boundary_cells
            if((coor_common_boundary_cell(1,ic)-xy(1))**2 + (coor_common_boundary_cell(2,ic)-xy(2))**2 < epsilon) then
              ! index_common_boundary(iface)= ic + icount_common_boundary_cell
              index_common_boundary(:,iface)= id_face(:,ic)
              find_common_boundary_cell(iface)= find_common_boundary_cell(iface) + 1
              ! common_boundary_xy(:,index_common_boundary(iface))= xy(:)
              icount= icount + 1
              write(*,*) 'processor ', MyProc, 'ipar',ipar,'index', index_common_boundary(:,iface)
              ! write(*,*) NR,index_common_boundary(iface),iface
              ! write(*,*) coor_common_boundary_cell(:,ic),xy(:)
              exit
            endif
          enddo
        endif
      enddo

      icount_common_boundary_cell=icount_common_boundary_cell+num_common_boundary_cells
      deallocate(coor_common_boundary_cell)
      deallocate(id_face)
    endif
  enddo



  if(maxval(find_common_boundary_cell)>1) then
    write(*,*) 'error in finding common boundary cells, one face has two common boundary neighbors, in processor', MyProc
    call MPI_STOP
  endif
  deallocate(find_common_boundary_cell)

  call MPI_AllReduce(icount,icount_sum, 1, MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(icount_sum /= icount_common_boundary_cell) then
    write(*,*) 'error: total number of common_boundary faces:', icount_common_boundary_cell,&
               'while the number of founded face:',icount_sum
    call MPI_STOP
  endif
  return
end subroutine