!!-----------------------------------------------------------
!! Mesh pre-process code, liwanai 2010.12.28, liwanai07@gmail.com
!! read in mesh file and partition region
!! output file for parallel computation
!!---------------------------------------------
    
    subroutine Partition
    use mainvar
    use ConnectFace
    implicit none
    integer,allocatable :: CellPar(:), &  !! NodePar(MaxNode), &
               ParCell(:,:), ParNode(:,:),NParC(:),NParN(:), NParC_Vir(:), &
               ParFace(:,:), NParF(:),NParF_Vir(:), &
               GloToLoc(:),CellGToL(:), &
               NBCell(:),BCell(:,:),icount(:)
    integer :: ncc(2),tarPar,tarCon,&
               NBod,ibod(10),bodF(200000,10),iType_Par(10),Peri(10),PeriLocToGlo(10)
    integer :: i,j,k,ic,icon,inbc,ipar,ipar1,ipar2,iv,nei,icell,iface,iface2,ipart, &
               nn,numvirc,ib,iL,inode,iperi,iver,np,ibelong,NL,NR,NListR,NListS, &
               NParN_Vir,NPeriodic,numface_Vir,iperi2,neic
    character*50 :: str,filename
    type(Connect), pointer :: tCon
    type(ConnectF),pointer :: tConF
    logical :: first
    logical,external :: contain
    INTEGER :: CSTAT, ESTAT
    CHARACTER(100) :: CMSG,command


    OPEN(87,FILE="par_info.in",STATUS='OLD')
    read(87,*) NPar
    close(87)
    write(*,*) 'Number of partitions:',Npar


    allocate(NConE(Npar));NConE=0
    allocate(ConE(MaxCon,Npar))
    do i=1,Npar
      do j=1,MaxCon
        ConE(j,i)%tarPar=0
        ConE(j,i)%tarCon=0
        ConE(j,i)%NListR=0
        ConE(j,i)%NListS=0
        ConE(j,i)%ListR(:,:)=0
        ConE(j,i)%ListS(:,:)=0
      enddo
    enddo

    allocate(NConF(Npar));NConF=0
    allocate(ConF(MaxCon,Npar))
    do i=1,Npar
      do j=1,MaxCon
        ConF(j,i)%tarPar=0
        ConF(j,i)%tarCon=0
        ConF(j,i)%NList=0
        ConF(j,i)%List(:)=0
      enddo
    enddo
    
    allocate(CellPar(NE));cellpar=0
    allocate(ParCell(NE,NPar));Parcell=0
    allocate(ParNode(NNT,NPar));ParNode=0
    allocate(NParC(NPar));NParC=0
    allocate(NParN(NPar));NparN=0
    allocate(NParC_Vir(NPar));NParC_Vir=0
    allocate(GloToLoc(NNT));GloToLoc=0
    allocate(CellGToL(NE));CellGToL=0
    allocate(NBCell(NPar));NBCell=0
    allocate(BCell(MaxVir,NPar));BCell=0
    allocate(icount(max(NE,NNT)));icount=0
    allocate(ParFace(NF,NPar));ParFace=0
    allocate(NParF(NPar));NParF=0
    allocate(NParF_Vir(NPar));NParF_Vir=0

  

    call metisgrid 

    !! read partition info
    if(NPar<10)then
      write(str,"(i1)") NPar
    elseif(NPar>=10 .and. NPar<100) then
      write(str,"(i2)") NPar
    elseif(NPar>=100 .and. NPar<1000) then
      write(str,"(i3)") NPar
    else
      write(*,*) 'Warning : Number of processors is no less than 1000!';stop
    endif

    command= trim("mpmetis grid_metis.in "//trim(str))
    
    CALL EXECUTE_COMMAND_LINE (command, EXITSTAT=ESTAT, &
           CMDSTAT=CSTAT, CMDMSG=CMSG)
    IF (CSTAT > 0) THEN
      write(*,*) "Command execution failed with error ", TRIM(CMSG)
    ELSE IF (CSTAT < 0) THEN
      write(*,*) "Command execution not supported"
    ELSE
      write(*,*) "Command completed with status", ESTAT
    END IF

    !! output to grid.in
    filename= trim("grid_metis.in.epart."//trim(str))
    
    OPEN(87,FILE=filename,STATUS='OLD')
    read(87,*) (CellPar(k),k=1,NE)
    close(87)
    CellPar(:)= CellPar+1
    write(*,*) 'finish reading partition info!'


    !! collect cells and nodes in every partition
    NParC(:)= 0
    NParN(:)= 0
    do i=1,NE
      
      !! add cell to this part
      ipart= CellPar(i)
      NParC(iPart)=NParC(ipart)+1
      ParCell( NParC(ipart),ipart )= i
      
      !! include 4 vertex of i, tetrahedron
      do iv=1, KFace(i)  ! 3
        !! N(1:3,i)
        nn= N(iv,i)
        do k=1,NparN(ipart)
        if( ParNode(k,ipart)== nn)exit
        enddo
        if( k<=NParN(ipart) )cycle

        !! add vertex to this part
        NparN(ipart)= NparN(ipart)+1
        ParNode(NparN(ipart),ipart)= nn
      enddo
    enddo

    !! determine the local index of the cell in the partion
    do ipar=1,NPar
      do ic=1,NParC(ipar)
        icell= ParCell(ic,ipar)
        CellGToL(icell)= ic
      enddo
    enddo
    
    write(*,*) 'finish collectting cells and nodes in every partition!'
    
    
!!----------------------------------------
!! find all boundary cells in each partition
!!----------------------------------------
  
    !! find the connectivities in all partitions
    NBCell(:)= 0
    do i=1, NF
      ncc(1:2)= Neighbor(1:2,i)
      if( ncc(2)<= 0 )  cycle !! boundary face
      if( CellPar(ncc(1))== CellPar(ncc(2)) )  cycle
      !! add face to Connected face ??

      !! add NL, NR to boundary cells, 会重复的??
      do ic=1, 2
        ipar= CellPar( ncc(ic) )

        !! check if ncc(ic) is chosen
        do iNbc=1,NBCell(ipar)
          if( BCell(iNbc,ipar)==ncc(ic) ) exit
        enddo
        if( iNbc<= NBCell(ipar) ) cycle

        NBCell(ipar)= NBCell(ipar)+1
        BCell(NBCell(ipar),ipar)= ncc(ic)
      enddo

      !! add this face to Connected face. repeated added problem? coressponding problem??
      !!  同时加的，应该不会造成不对应吧? 最终可以检查下
      do ic=1,2
        ipar= CellPar( ncc(ic) )
        tarPar= CellPar( ncc( mod(ic,2)+1 ) )

        !! which ConF(icon,ipar) should I add to
        do icon=1,NConF(ipar)
          if( ConF(icon,ipar)%tarPar== tarPar) then
            tConF=> ConF(icon,ipar)
            tConF%NList = tConF%NList+ 1
            tConF%List( tConF%NList )= i
            exit
          endif
        enddo
        
        !! create new connected face in this partion
        if(icon>NConF(ipar))then
          NConF(ipar)= NConF(ipar)+ 1;  icon= NConF(ipar)
          tConF=> ConF(icon,ipar)
          tConF%tarPar= tarPar
          ! tConF%tarCon= ?
          tConF%NList= 1
          tConF%List(1)= i
        endif
      enddo

    enddo


    !! tell the Connected face in every partition, tarCon
    do ipar=1, NPar
    do icon=1, NConF(ipar)
      tarPar= ConF(icon,ipar)%tarPar

      do tarCon=1, NConF(tarPar)
      if( ConF(tarCon,tarPar)%tarPar== ipar )exit
      enddo
      if( tarCon>NConF(tarPar) )then
        write(*,*) 'error in tarCon in connected face'; stop
      endif

      !! add tarCon
      ConF(icon,ipar)%tarCon= tarCon
    enddo
    enddo


    !! sorts BCell's virtual neighbors
    do ipar=1, NPar
      numVirC= NParC(ipar)
    do ic=1,NBCell(ipar)
      icell= BCell(ic,ipar)
      do nei=1, Neigh(icell)
        neic= NeighC(nei,icell)
        tarpar= CellPar(neic)
        if( tarpar == ipar )cycle

        !! check if neic has virtual cells? 这个还是没有做好???
        !! add and sort neic
        do icon=1, NConE(ipar)
        if( ConE(icon,ipar)%TarPar == tarpar )then
          do iL=1,ConE(icon,ipar)%NListR
          if( ConE(icon,ipar)%ListR(0,iL)== neic ) goto 111
          enddo
        endif
        enddo
111    if( icon<= NConE(ipar) )   cycle

        
        !! add and sort neic
        do icon=1, NConE(ipar)
        if( ConE(icon,ipar)%TarPar == tarpar )then
          ! add neic to this connect, which info needed ?
          NListR= ConE(icon,ipar)%NListR +1; numVirC= numVirC+1
          ConE(icon,ipar)%NListR= NListR
          ConE(icon,ipar)%ListR(0:2,NListR)= (/neic, numVirC, CellGToL(neic) /)

          !! add virtual cell to ParCell(:,:)
          ParCell(numVirC,ipar)= neic !! give global index to ParCell
          exit
        endif
        enddo
        if(icon<=NConE(ipar)) cycle

        !! creat new connect and put neic to this connect
        NConE(ipar)= NConE(ipar)+1;  numVirC= numVirC+1
        ConE(NConE(ipar),ipar)%TarPar= tarpar
        ConE(NConE(ipar),ipar)%NListR= 1
        ConE(NConE(ipar),ipar)%ListR(0:2,1)= (/ neic, numVirC, CellGToL(neic) /)
                                            !! (/ global index, local index in this par, local index in target par /)
        !! add virtual cell to ParCell(:,:)
        ParCell(numVirC,ipar)= neic !! give global index to ParCell
      enddo !! --- nei
    enddo

      NParC_Vir(ipar)= numVirC
      write(*,*) 'finish adding vitrual cells and nodes in partion:',ipar
    enddo


    !! for each connect, tell the target parition my required cells
    do ipar=1, NPar
    do icon=1, NConE(ipar)
      tarPar= ConE(icon,ipar)%tarPar
      !! tarCon?
      do tarCon=1,NConE(tarPar)
      if( ConE(tarCon,tarPar)%TarPar == ipar )exit
      enddo
      if( tarCon>NConE(tarPar) )then
        write(*,*) "error in tarCon!"; stop
      endif
      
      ConE(icon,ipar)%tarCon= tarCon
      !! tell tarCon at tarPar my requirement
      NListS= ConE(icon,ipar)%NListR
      ConE(tarCon,tarPar)%NListS= NListS
      ConE(tarCon,tarPar)%ListS(1,1:NListS)= ConE(icon,ipar)%ListR(2,1:NListS)
      ConE(tarCon,tarPar)%ListS(2,1:NListS)= ConE(icon,ipar)%ListR(1,1:NListS)

    enddo
    enddo

!!=================================================
!! compute the number of faces in every partion
!!-------------------------------------------------

    NParF=0
    do i=1,NF
      NL=neighbor(1,i)
      NR=neighbor(2,i)
      ipar1=cellpar(NL)
      NParF(ipar1)=NParF(ipar1)+1
      ParFace(NParF(ipar1),ipar1)=i
      if(NR<0) cycle
      if(Fproperty(i)==100) cycle
      ipar2=cellpar(NR)
      if(ipar1 /= ipar2) then 
        NParF(ipar2)=NParF(ipar2)+1
        ParFace(NParF(ipar2),ipar2)=i
      endif
    enddo
    
    do ipar=1,Npar
        numface_Vir=NParF(ipar)
        first= .true.
    do ic=NParC(ipar),NParC_Vir(ipar)
      icell= ParCell( ic,ipar )
      do iface=1,KFace(icell)
         iface2=NeighF(iface,icell)
         NL=neighbor(1,iface2)
         NR=neighbor(2,iface2)
         if(NR>0) then
           if(Fproperty(iface2)==100) then
            if(cellpar(NL)==ipar) cycle
           else
            if(cellpar(NL)==ipar .or. cellpar(NR)==ipar) cycle
           endif
         else
           if(cellpar(NL)==ipar) cycle  
         endif
      
      if(first) then  
         numface_Vir=numface_Vir+1
         ParFace(numface_Vir,ipar)=iface2
         first= .false.
      else
        if( .not. contain( iface2, ParFace( (NParF(ipar)+1):numface_Vir,ipar), numface_Vir-NParF(ipar) ) ) then
          numface_Vir=numface_Vir+1
          ParFace(numface_Vir,ipar)=iface2
        endif
      endif
      enddo
    enddo
      NParF_Vir(ipar)=numface_Vir
    enddo
     
    
    !do ipar=1,Npar
    !    numface_Vir=NParF(ipar)
    !do ic=NParC(ipar),NParC_Vir(ipar)
    !  icell= ParCell( ic,ipar )
    !  do iface=1,KCFace(icell)
    !    if( .not. contain( NeighF(iface,icell), ParFace(1:numface_Vir,ipar), numface_Vir ) ) then
    !      numface_Vir=numface_Vir+1
    !      ParFace(numface_Vir,ipar)=NeighF(iface,icell)
    !    endif
    !  enddo
    !enddo
    !  NParF_Vir(ipar)=numface_Vir
    !enddo

!!=================================================
!! output every part, grid.in
!!-----------------------------------------
    do ipar=1,NPar

    !! 本区包含的点的局部编号, global->local
    GloToLoc(:)=-1
    do inode=1,NParN(ipar)  !! Local: in, Global: nn
      nn= ParNode(inode,ipar)
      GloToLoc(nn)= inode
    enddo
    
    !! add virtual vertex
    NParN_Vir= NParN(ipar)
    do icon=1, NConE(ipar)
    do iL=1, ConE(icon,ipar)%NListR
      ic= ConE(icon,ipar)%ListR(0,iL)
      do iver=1,KFace(ic)  !! 3
        nn= N(iver,ic)
        !! nn是否在ParNode(:,ipar)? need to cycle all vertex in this? No
        if( GloToLoc(nn)>0 )cycle
        !! add nn to this partition
        NParN_Vir= NParN_Vir+1
        ParNode(NParN_Vir,ipar)= nn
        GloToLoc(nn)= NParN_Vir
      enddo
    enddo
    enddo
    
     !! find boundary including the periodic boundary
    ibod(:)=0; bodF(:,:)=0; NBod=0; Nperiodic=0
    do ib=1, NBoundary
      first= .true.
    do np=1, ibound(ib)
      iface= BoundF(np,ib)
          ibelong=0
          !if(CellPar(Neighbor(1,iface))==ipar) ibelong=1
          !if(Neighbor(2,iface)>0) then
          !   if(FProperty(iface)== 100 .and. CellPar(Neighbor(2,iface))==ipar) ibelong=1
          !endif
          
          !! inner faces
          if(CellPar(Neighbor(1,iface))==ipar) ibelong=1
          !! virtual faces
          if( contain( iface, ParFace( (NParF(ipar)+1):NParF_Vir(ipar),ipar), NParF_Vir(ipar)-NParF(ipar) ) ) ibelong=1
          
      if(ibelong>0) then
      !if( CellPar(Neighbor(1,iface))==ipar &  !! iface belong to this parition
      !    .or. ( FProperty(iface)== 100 .and. CellPar(Neighbor(2,iface))==ipar ))then  !! periodic neighbor in this partion
        if( first ) then
          NBod= NBod+1  !! 有问题????
          iType_Par(NBod)= iType(ib) !!ib, modified by Qian Wang
          if(iType_Par(NBod)==100) then
            Nperiodic=Nperiodic+1
            Peri(Nperiodic)=NBod
            PeriLocToGlo(Nperiodic)=ib
          endif
          first= .false.
        endif
        ibod(NBod)= ibod(NBod)+1
        bodF(ibod(NBod),NBod)= iface
      endif
    enddo
    enddo
    

    if(ipar<10)then
      write(str,'(i1)') ipar
    elseif(ipar>=10 .and. ipar<100) then
      write(str,'(i2)') ipar
    elseif(ipar>=100 .and. ipar<1000) then
      write(str,'(i3)') ipar
    else
      write(*,*) 'Warning : Number of processors is no less than 1000!';stop
    endif
    

    !! output to grid.in
    filename= trim("grid"//trim(str)//".in")
    open(unit=22,file=filename)
    write(22,*) NParC(ipar),NParC_Vir(ipar),NParN(ipar),NParN_Vir,NParF(ipar),NParF_Vir(ipar)
    do ic=1,NParC_Vir(ipar)  !! NParC(ipar)   !! Local: ic, Global: icell
      icell= ParCell( ic,ipar )
      write(22,*) KFace(icell),( GloToLoc( N(k,icell) ), k=1,KFace(icell) )
    enddo

    do inode=1,NParN_Vir  !! NParN(ipar)
      nn= ParNode(inode,ipar)
      write(22,*) (Coor(k,nn),k=1,2), (IsCommenBoundaryNode(k,nn),k=1,2)
    enddo


    !! output each connect arrays: virtual cells and target cells
    write(22,*) NCone(ipar),  ',        **** connected cells'
    
    do icon=1, NConE(ipar)
      tCon => ConE(icon,ipar)
      write(22,*) tCon%tarPar,tCon%tarCon, tCon%NListR,tCon%NListS
      do iL=1, tCon%NListR
        write(22,*) (tCon%ListR(k,iL),k=1,2)
      enddo
      do iL=1, tCon%NListS
        write(22,*) (tCon%ListS(k,iL),k=1,2)
      enddo
    enddo

    !!! output connected faces
    !write(22,*) NconF(ipar), ',         *** connected face'
    !do icon=1, NConF(ipar)
    !  tConF => ConF(icon,ipar)
    !  write(22,*) tConF%tarPar,tConF%tarCon,tConF%NList, ' , '
    !  do iL=1,tConF%NList
    !  !! output vertex local index of the face
    !    write(22,*) KFNode(tConF%List(iL)),(GloToLoc( FNode(k,tConF%List(iL)) ), k=1,KFNode(tConF%List(iL)) )
    !  enddo
    !enddo
    
    
    !! output boundary
    write(22,*) NBod,maxval(ibod(:))
    do ib=1, NBod
      write(22,*) iType_par(ib), ibod(ib), ',    *'
      do np=1,ibod(ib)
        iface= BodF(np,ib)
        write(22,*) ( GloToLoc(FNode(k,iface)), k=1,2)
      enddo
    enddo
    
    
!! write periodic boudanry info

      write(22,*) Nperiodic
      if(Nperiodic>0) then
      do iperi=1,Nperiodic
          
          do icon=1,Ncon
              if(PeriLocToGlo(iperi)==Con(1,icon)) exit
          enddo
          if(icon>Ncon) then
              write(*,*) 'error in finding the periodic boudanry in Partion:',ipar
              stop
          endif
          
          do iperi2=1,Nperiodic
          if(PeriLocToGlo(iperi2)==Con(2,icon)) exit
          enddo
          if(iperi2>Nperiodic) then
              write(*,*) 'error in periodic boudanry matching in Partion:',ipar
              stop
          endif

          write(22,*) Peri(iperi),Peri(iperi2), (dxyz(k,icon),k=1,2)
      enddo
      endif

    close(22)
    
    
    !! output to tecplot

    ! filename= trim("grid"//trim(str)//".dat")
    ! open(unit=22,file=filename)
    ! write(22,*) 'variables="x","y","ipar"'
    ! write(22,*) 'ZONE N=',NParN(ipar),',E=',NParC(ipar),', F=FEPOINT, ET=QUADRILATERAL' 
    ! do inode=1,NParN(ipar)
    !   nn= ParNode(inode,ipar)
    !   write(22,*) (Coor(k,nn),k=1,2),ipar
    ! enddo
    
    ! do ic=1,NParC(ipar)   !! Local: ic, Global: icell
    !   icell= ParCell( ic,ipar )
    !   write(22,*) (GloToLoc( N(k,icell)),k=1,4)
    ! enddo
    ! close(22)

    ! !! output to tecplot with virtual cells
    ! filename= trim("grid_Vir"//trim(str)//".dat")
    ! open(unit=22,file=filename)
    ! write(22,*) 'variables="x","y","ipar"'
    ! write(22,*) 'ZONE N=',NParN_Vir,',E=',NParC_Vir(ipar),', F=FEPOINT, ET=QUADRILATERAL'
    ! do inode=1,NParN_Vir
    !   nn= ParNode(inode,ipar)
    !   write(22,*) (Coor(k,nn),k=1,2),ipar
    ! enddo
    
    ! do ic=1,NParC_Vir(ipar)   !! Local: ic, Global: icell
    !   icell= ParCell( ic,ipar )
    !   write(22,*) (GloToLoc( N(k,icell)),k=1,4)
    ! enddo
    ! close(22)

    enddo
    
    
    !! check send and recieve vortual grids
    
    if(sum(NParC(:)) /= NE) then
      write(*,*) 'error in the total number of the cells in the partions '
      stop
    endif
    
    !! avoid repeatation of the cells
    icount(:)=0
    do ipar=1,Npar
    do ic=1,NParC(ipar)   !! Local: ic, Global: icell
      icell= ParCell( ic,ipar )
      icount(icell)=icount(icell)+1
    enddo
    enddo
    
    do i=1,NE
    if(icount(i) /= 1) then
      write(*,*) 'error in cell distribution '
      stop
    endif
    enddo
    
    !! avoid repeatation of the nodes
    icount(:)=0
    do ipar=1,Npar
    do inode=1,NParN(ipar)
      nn= ParNode(inode,ipar)
      icount(nn)=icount(nn)+1
    enddo
    enddo
    
    do ipar=1,NPar
    do icon=1,NConF(ipar)      
      tConF=> ConF(icon,ipar)
    do iL=1,tConF%NList
      iface=tConF%List( iL )
      if(Fproperty(iface) /= 100 ) then
      do inode=1,2
        icount(Fnode(inode,iface))=1
      enddo
      endif
    enddo
    enddo
    enddo
    
    do i=1,NNT
    if(icount(i) /= 1) then
      write(*,*) 'error in node distribution '
      stop
    endif
    enddo
    
    !! check the send arraries
    icount=0
    
    do i=1,NE
     do nei=1,Neigh(i)
     if(CellPar(NeighC(nei,i))==CellPar(i))  icount(i)=icount(i)+1
     enddo
    enddo
    
    do ipar=1,NPar
    do icon=1, NConE(ipar)
      tCon=> ConE(icon,ipar)
      do ic=1,tCon%NListS
        icell=Parcell(tCon%ListS(1,ic),ipar)
        do nei=1,Neigh(icell)
          if(CellPar(NeighC(nei,icell))==tCon%tarPar)  icount(icell)=icount(icell)+1
        enddo
      enddo
    enddo
    enddo
    
    do i=1,NE
    if(icount(i) /= Neigh(i)) then
      write(*,*) 'error in send neighbor '
      stop
    endif
    enddo
    
    !! check recieve arraries
    icount=0
    
    do i=1,NE
     do nei=1,Neigh(i)
     if(CellPar(NeighC(nei,i))==CellPar(i))  icount(i)=icount(i)+1
     enddo
    enddo
    
    do ipar=1,NPar
    do icon=1, NConE(ipar)
      tCon=> ConE(icon,ipar)
      do ic=1, tCon%NListR
        icell= tCon%ListR(0,ic)
        do nei=1,Neigh(icell)
        if( CellPar( NeighC(nei,icell) )==ipar ) icount( NeighC(nei,icell) )=icount( NeighC(nei,icell) )+1
        enddo
      enddo
    enddo
    enddo
    
    do i=1,NE
    if(icount(i) /= Neigh(i)) then
      write(*,*) 'error in recieve neighbor '
      stop
    endif
    enddo
    
    write(*,*) 'success !!!'
    
    
    deallocate(CellPar)
    deallocate(ParCell)
    deallocate(ParNode)
    deallocate(NParC)
    deallocate(NParN)
    deallocate(NParC_Vir)
    deallocate(GloToLoc)
    deallocate(CellGToL)
    deallocate(NBCell)
    deallocate(BCell)
    deallocate(icount)
    deallocate(ParFace)
    deallocate(NParF)
    deallocate(NParF_Vir)

    deallocate(NConE)
    deallocate(ConE)
    deallocate(NConF)
    deallocate(ConF)
    
    return
    end subroutine
