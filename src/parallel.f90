!!***************************************************************
!! data transformation subroutines for parallel computing
!!***************************************************************


module 	Parallel_Var
    implicit none

    !! parallel info
    integer,save :: NPar
    integer,save :: max_NListR,max_NListS
    integer,save :: myProc
    character*50 :: strProc,filename
    integer,parameter :: bufferSize=20000000
    real*8 :: bufferForISend(bufferSize)
    !! connectivity info
    type Connect
      integer :: tarPar,tarCon, NListR, NListS   !! NListR==NListS ??
      integer,allocatable :: ListR(:,:),ListS(:,:) !! ListR(2,NListR),ListS(2,NListS)
       !! List(2,NList): 1: this partition local index(virtual); 2: target local index
    end type Connect
    
    integer,save :: NConE
    type(Connect),allocatable,target :: ConE(:)
end module

!!----------------------------------------------
!! init parallelization
!!-------------------------------
subroutine InitParallel
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer  myrank,ierr, numprocess
    character*50,external:: getstring

    call MPI_INIT( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr)
    MyProc= myrank+1
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocess, ierr)
    print *, "Process ", myrank, " of ", numprocess, " is alive"
    if( numProcess/= NPar )then
      write(*,*) 'error, parallel setting file /= run parameter'
      stop
    endif

    !! allocate isend buffer size
    call MPI_Buffer_Attach(bufferForISend,bufferSize, ierr)

    strProc= getstring(MyProc)

    return
end subroutine


subroutine MPI_STOP
  implicit none
  include 'mpif.h'
  integer ierr
    stop
    ! call MPI_Finalize(ierr)
end subroutine



!!-----------------------------------------------
!! physical data transform on virtual grids
!!-----------------------------------------
subroutine Connect_CellAverage
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer :: icon,ic,icell,  TargetProc,tag,numGrid,NListS
    integer :: sendq,STATUS(MPI_STATUS_SIZE),ierr
    real*8,allocatable :: bufs(:,:),bufr(:,:)
    type(Connect),pointer :: tCon

    !! improve efficiency of packing, send data
    !! send data
    do icon=1, NConE
      tCon=> ConE(icon)
      TargetProc= tCon%tarPar -1 !!?????? 有问题 acutal process index
      tag       = tCon%tarCon
      numGrid   = NVar* tCon%NListS
      NListS= tCon%NListS
      allocate(bufs(NVar,NListS))
      do ic=1,NListS
        bufs(1:NVar,ic)= PA(1:NVar,tCon%ListS(1,ic))
      enddo
      call MPI_IBSend(bufs ,numGrid,MPI_DOUBLE_PRECISION,TargetProc,tag, &
                      MPI_COMM_WORLD, sendq, ierr)
      call MPI_Wait(sendq,status,ierr)

      deallocate(bufs)
    enddo

    do icon=1, NConE
      tCon=> ConE(icon)
      targetProc= tCon%tarPar -1 !! ???? 有问题 actual process index
      tag       = icon  !! tCon%tarCon
      numGrid   = NVar*tCon%NListR
      allocate(bufr(NVar,tCon%NListR) )
      call MPI_Recv(bufr,numGrid,MPI_DOUBLE_PRECISION,targetProc,tag, &
                 MPI_COMM_WORLD,status,ierr)

      !! set virtual grids' value
      do ic=1, tCon%NListR
        icell= tCon%ListR(1,ic)
        PA(1:NVar,icell)= bufr(1:NVar,ic)
      enddo

      deallocate(bufr)
    enddo

!    write(*,*) myProc,'out connect'

end subroutine


!!-----------------------------------------------
!! physical data transform on virtual grids
!!-----------------------------------------
subroutine Connect_LUSGS
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer :: icon,ic,icell,  TargetProc,tag,numGrid,NListS
    integer :: sendq,STATUS(MPI_STATUS_SIZE),ierr
    real*8,allocatable :: bufs(:,:),bufr(:,:)
    type(Connect),pointer :: tCon


    !! improve efficiency of packing, send data
    !! send data
    do icon=1, NConE
      tCon=> ConE(icon)
      TargetProc= tCon%tarPar -1 !!?????? 有问题 acutal process index
      tag       = tCon%tarCon
      numGrid   = NVar* tCon%NListS
      NListS= tCon%NListS
      allocate(bufs(NVar,NListS))
      do ic=1,NListS
        bufs(1:NVar,ic)= hWA(1:NVar,tCon%ListS(1,ic))
      enddo
      call MPI_IBSend(bufs ,numGrid,MPI_DOUBLE_PRECISION,TargetProc,tag, &
                      MPI_COMM_WORLD, sendq, ierr)
      call MPI_Wait(sendq,status,ierr)

      deallocate(bufs)
    enddo

    do icon=1, NConE
      tCon=> ConE(icon)
      targetProc= tCon%tarPar -1 !! ???? 有问题 actual process index
      tag       = icon  !! tCon%tarCon
      numGrid   = NVar*tCon%NListR
      allocate(bufr(NVar,tCon%NListR) )
      call MPI_Recv(bufr,numGrid,MPI_DOUBLE_PRECISION,targetProc,tag, &
                 MPI_COMM_WORLD,status,ierr)

      !! set virtual grids' value
      do ic=1, tCon%NListR
        icell= tCon%ListR(1,ic)
        hWA(1:NVar,icell)= bufr(1:NVar,ic)
      enddo

      deallocate(bufr)
    enddo

!    write(*,*) myProc,'out connect'

end subroutine

!!--------------------------------------------------
!! necessary when have limiting
!!--------------------------------------------------
    subroutine Connect_Gradient
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer :: icon,ic,icell,k,  TargetProc,tag,numGrid,NListS
    integer :: sendq,STATUS(MPI_STATUS_SIZE),ierr
    real*8,allocatable :: bufs(:,:,:),bufr(:,:,:)
    type(Connect),pointer :: tCon

    !! improve efficiency of packing, send data
    !! send data
    do icon=1, NConE
      tCon=> ConE(icon)
      TargetProc= tCon%tarPar -1 !!?????? 有问题 acutal process index
      tag       = tCon%tarCon
      numGrid   = MaxOrder*NVar*tCon%NListS
      NListS= tCon%NListS
      allocate(bufs(MaxOrder,NVar,NListS))
      do ic=1, NListS
      do k =1, NVar  !! 5   !4
        bufs(1:MaxOrder,k,ic)= gradC(1:MaxOrder,k,tCon%ListS(1,ic))
      enddo
      enddo
      call MPI_IBSend(bufs ,numGrid,MPI_DOUBLE_PRECISION,TargetProc,tag, &
                      MPI_COMM_WORLD, sendq, ierr)
      call MPI_Wait(sendq,status,ierr)

      deallocate(bufs)
    enddo

    do icon=1, NConE
      tCon=> ConE(icon)
      targetProc= tCon%tarPar -1 !! ???? 有问题 actual process index
      tag       = icon  !! tCon%tarCon
      numGrid   = MaxOrder*NVar* tCon%NListR
      allocate( bufr(MaxOrder,NVar,tCon%NListR) )
      call MPI_Recv(bufr,numGrid,MPI_DOUBLE_PRECISION,targetProc,tag, &
                    MPI_COMM_WORLD,status,ierr)

      !! set virtual grids' value
      do ic=1, tCon%NListR
        icell= tCon%ListR(1,ic)
        do k=1, NVar ! 5 !4
        gradC(1:MaxOrder,k,icell)= bufr(1:MaxOrder,k,ic)
        enddo
      enddo

      deallocate(bufr)
    enddo

    end subroutine

!!--------------------------------------------------
!! necessary when using H-Entropy Correction
!!--------------------------------------------------
    subroutine Connect_Deltaeig
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer :: icon,ic,icell,k,  TargetProc,tag,numGrid,NListS
    integer :: sendq,STATUS(MPI_STATUS_SIZE),ierr
    real*8,allocatable :: bufs(:),bufr(:)
    type(Connect),pointer :: tCon

    !! improve efficiency of packing, send data
    !! send data
    do icon=1, NConE
      tCon=> ConE(icon)
      TargetProc= tCon%tarPar -1 !!?????? 有问题 acutal process index
      tag       = tCon%tarCon
      numGrid   = tCon%NListS
      NListS= tCon%NListS
      allocate(bufs(NListS))
      do ic=1, NListS
        bufs(ic)= deltaeig(tCon%ListS(1,ic))
      enddo
      call MPI_IBSend(bufs ,numGrid,MPI_DOUBLE_PRECISION,TargetProc,tag, &
                      MPI_COMM_WORLD, sendq, ierr)
      call MPI_Wait(sendq,status,ierr)

      deallocate(bufs)
    enddo

    do icon=1, NConE
      tCon=> ConE(icon)
      targetProc= tCon%tarPar -1 !! ???? 有问题 actual process index
      tag       = icon  !! tCon%tarCon
      numGrid   = tCon%NListR
      allocate( bufr(tCon%NListR) )
      call MPI_Recv(bufr,numGrid,MPI_DOUBLE_PRECISION,targetProc,tag, &
                    MPI_COMM_WORLD,status,ierr)

      !! set virtual grids' value
      do ic=1, tCon%NListR
        icell= tCon%ListR(1,ic)
        deltaeig(icell)= bufr(ic)
      enddo

      deallocate(bufr)
    enddo

    end subroutine


!!---------------------------------------
!!  checking connectivity coordinate
!!---------------------------------------

subroutine Connect_CellCenter
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer :: icon,ic,icell,  TargetProc,tag,numGrid,NListS,k
    integer :: sendq,STATUS(MPI_STATUS_SIZE),ierr
    real*8,allocatable :: bufs(:,:),bufr(:,:)
    type(Connect),pointer :: tCon

    !! improve efficiency of packing, send data
    !! send data
    do icon=1, NConE
      tCon=> ConE(icon)
      TargetProc= tCon%tarPar -1 !!?????? 有问题 acutal process index
      tag       = tCon%tarCon
      numGrid   = 2*tCon%NListS
      NListS= tCon%NListS
      allocate(bufs(2,NListS))
      do ic=1,NListS
        bufs(1:2,ic)= CellXY(1:2,tCon%ListS(1,ic))
      enddo
      call MPI_IBSend(bufs ,numGrid,MPI_DOUBLE_PRECISION,TargetProc,tag, &
                      MPI_COMM_WORLD, sendq, ierr)
      call MPI_Wait(sendq,status,ierr)

      deallocate(bufs)
    enddo

    do icon=1, NConE
      tCon=> ConE(icon)
      targetProc= tCon%tarPar -1 !! ???? 有问题 actual process index
      tag       = icon  !! tCon%tarCon
      numGrid   = 2*tCon%NListR
      allocate(bufr(2,tCon%NListR) )
      call MPI_Recv(bufr,numGrid,MPI_DOUBLE_PRECISION,targetProc,tag, &
                 MPI_COMM_WORLD,status,ierr)

      !! compare the center coordinate
      do ic=1, tCon%NListR
        icell= tCon%ListR(1,ic)
        if(icell<=NE) then
            write(*,*) 'virtual index lower than inner cells in process:',MyProc 
            stop
        endif
!        write(*,*) myProc,bufr(1,ic),CellXY(1,icell)-bufr(1,ic),CellXY(2,icell)-bufr(2,ic)
        if( abs(CellXY(1,icell)-bufr(1,ic)) + abs(CellXY(2,icell)-bufr(2,ic)) > 1.e-8)then
          !if( abs( sum(abs(CellXY(:,icell)-bufr(:,ic)))-0.2) >1.e-4 )then
          write(*,*) 'error in coordinate in',myProc
          write(*,*) 'the error is located at icon:',icon
          write(*,*) 'target partion is ', targetProc+1
          write(*,*) 'the error of the coordinates are',(CellXY(k,icell)-bufr(k,ic),k=1,2)
          stop
          !endif
        endif
      enddo

      deallocate(bufr)
    enddo

end subroutine

