subroutine metisgrid
  use mainvar
  use ConnectFace
  implicit none
  integer :: igridtype,i,k,np,ib,icon
  
  !OPEN(5,FILE='grid_db_240.in') 
  !write(5,'(2I)')  NE,NNT
  !
  !do i=1,NE
  !write(5,*) Kface(i),( N(k,i), k=1,KFace(i) )
  !enddo
  !
  !do np=1,NNT
  !write(5,*) COOR(1,NP),COOR(2,NP)
  !enddo
  !
  !  write(5,*) NBoundary,maxval(ibound(:)),  '    **boundaryinfo'
  !  
  !  do ib=1, NBoundary
	 ! write(5,*) itype(ib),ibound(ib), ' *'
  !    do np=1,ibound(ib)
	 ! write(5,*) BoundN(1,np,ib),BoundN(2,np,ib)
	 ! enddo
  !  enddo
  !  
  !  write(5,*) NCon
  !  do icon=1,NCon
  !    write(5,*) (Con(k,icon),k=1,2), (dxyz(k,icon),k=1,2)
  !  enddo
  !
  !close(5)
  if(MaxNei==4) then
  igridtype=4 !! 1:tri,2:tetra,3:hexa,4:quad
  else
  igridtype=1 !! 1:tri,2:tetra,3:hexa,4:quad
  endif
  OPEN(111,FILE='grid_metis.in') 
  write(111,*)  NE !,igridtype
  do i=1,NE
  write(111,*) ( N(k,i), k=1,KFace(i) )
  enddo
  close(111)
end subroutine