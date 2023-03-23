
!!==============================
!!  output : cell-center format
!!------------------------------
SUBROUTINE output_vtk(jout)
    use mainvar
    use Parallel_Var
    implicit none
    double precision,allocatable:: RHO(:),U(:),V(:),PRE(:),MA(:),TEM(:),Y(:,:)
    double precision:: Cp, gama, Rcpcv, ac        
    integer :: i,jout,k, Nsize,ks
    integer, parameter :: IOfile= 121
    character*50  str,str2
    character*150 str3
    character*50,external:: getstring
    
    allocate(RHO(NE));RHO=0.
    allocate(U(NE));U=0.
    allocate(V(NE));V=0.
    allocate(PRE(NE));PRE=0.
    allocate(MA(NE));MA=0.
    allocate(TEM(NE));TEM=0.
    allocate(Y(Nvar-Nflow+1,NE));Y=0.
    ! allocate(YSUM(NE));YSUM=0.
    ! allocate(PROD_SUM(NE));PROD_SUM=0.

    do i=1,NE
      pre(i)= PA(1,i)
      u(i)  = PA(2,i)
      v(i)  = PA(3,i)
      tem(i)= PA(4,i)
      Y(1:Nvar-Nflow,i)= PA(nflow+1:NVar,i)

      Y(Nvar-Nflow+1,i)= 1. - sum(Y(1:Nvar-Nflow,i))
      
      call ComputeGasParameter(Y(:,i), tem(i), Cp, gama, Rcpcv)
      rho(i)= pre(i)/(Rcpcv*tem(i))
      ac= sqrt(gama*Rcpcv*tem(i))
      Ma(i)=  sqrt(u(i)**2 + v(i)**2)/ac

      ! YSUM(i)= sum(Y(:,i))
      ! PROD_SUM(i)= sum(Production_rate(:,i))

      if((pre(i)<0 .or. rho(i)<0) .or. tem(i)<0) then
        write(*,*) 'Error: negative density/pressure/temperature:', rho(i), pre(i), tem(i)
        write(*,*) 'Cell index:',i
        write(*,*) 'Location of the cell:',CellXY(1,i),CellXY(2,i)
        call MPI_STOP
      endif

      if(abs(Production_rate(1,i))<1.0e-15) Production_rate(:,i)=0.

      do k=1,Nvar-Nflow
        if(abs(Y(k,i))<1.0e-15) Y(k,i)=0.
      enddo  
    enddo

      !!  create a .vtk file
      filename= trim("output"//trim(getstring(jout))//"_"//trim(strProc)//".vtk")
      OPEN(IOfile,FILE=filename) 

      !! header
      str=trim("# vtk DataFile Version 1.0")
      write(IOfile,'(a50)') str
      str=trim("Unstructured Grid Format")
      write(IOfile,'(a50)') str
      str=trim("ASCII")
      write(IOfile,'(a50)') str
      
      str=trim("DATASET UNSTRUCTURED_GRID")
      write(IOfile,'(a50)') str

      !! coordinates
      str= trim("POINTS "//trim(getstring(NNT))//" double")
      write(IOfile,'(a50)') str
      do i=1,NNT
        write(str3,'(3e15.8)') Coor(1,i), Coor(2,i), 0.
        write(IOfile,'(a)') adjustl(trim(str3))
      enddo

      !! Cell nodes
      Nsize=sum(kface(1:NE)) + NE
      str= trim("CELLS "//trim(getstring(NE))//" "//trim(getstring(Nsize)))
      write(IOfile,'(a50)') str
      do i=1,NE
        if (kface(i)==3) then
          write(str3,'(4I10)') kface(i), (N(k,i)-1,k=1,kface(i))
          write(IOfile,'(a)') adjustl(trim(str3))
        else
          write(str3,'(5I10)') kface(i), (N(k,i)-1,k=1,kface(i))
          write(IOfile,'(a)') adjustl(trim(str3))
        endif
        
      enddo

      !! Cell types: tri:5, quad:9
      str= trim("CELL_TYPES "//trim(getstring(NE)))
      write(IOfile,'(a50)') str
      do i=1,NE
        if (kface(i)==3) then
          write(str3,'(I1)') 5
          write(IOfile,'(a)') adjustl(trim(str3))
        else
          write(str3,'(I1)') 9
          write(IOfile,'(a)') adjustl(trim(str3))
        endif
      enddo

      !! cell-center data
      str= trim("CELL_DATA "//trim(getstring(NE)))
      write(IOfile,'(a50)') str

      str= trim("SCALARS rho double")
      write(IOfile,'(a50)') str
      str=trim("LOOKUP_TABLE default")
      write(IOfile,'(a50)') str
      do i=1,NE
        write(str3,'(e20.8)') rho(i)
        write(IOfile,'(a)') adjustl(trim(str3))
      enddo

      str= trim("SCALARS u double")
      write(IOfile,'(a50)') str
      str=trim("LOOKUP_TABLE default")
      write(IOfile,'(a50)') str
      do i=1,NE
        write(str3,'(e20.8)') U(i)
        write(IOfile,'(a)') adjustl(trim(str3))
      enddo

      str= trim("SCALARS v double")
      write(IOfile,'(a50)') str
      str=trim("LOOKUP_TABLE default")
      write(IOfile,'(a50)') str
      do i=1,NE
        write(str3,'(e20.8)') V(i)
        write(IOfile,'(a)') adjustl(trim(str3))
      enddo

      str= trim("SCALARS p double")
      write(IOfile,'(a50)') str
      str=trim("LOOKUP_TABLE default")
      write(IOfile,'(a50)') str
      do i=1,NE
        write(str3,'(e20.8)') pre(i)
        write(IOfile,'(a)') adjustl(trim(str3))
      enddo

      str= trim("SCALARS T double")
      write(IOfile,'(a50)') str
      str=trim("LOOKUP_TABLE default")
      write(IOfile,'(a50)') str
      do i=1,NE
        write(str3,'(e20.8)') tem(i)
        write(IOfile,'(a)') adjustl(trim(str3))
      enddo

      str= trim("SCALARS Ma double")
      write(IOfile,'(a50)') str
      str=trim("LOOKUP_TABLE default")
      write(IOfile,'(a50)') str
      do i=1,NE
        write(str3,'(e20.8)') Ma(i)
        write(IOfile,'(a)') adjustl(trim(str3))
      enddo

      do ks=1, Nvar-Nflow+1
        ! write(IOfile,*) 
        str=trim("SCALARS Y"//trim(getstring(ks))//" double")
        write(IOfile,'(a)') adjustl(trim(str))
        str=trim("LOOKUP_TABLE default")
        write(IOfile,'(a)') adjustl(trim(str))
        do i=1,NE
          write(str3,'(es30.12)') Y(ks,i)
          write(IOfile,'(a)') adjustl(trim(str3))
        enddo
      enddo

      ! str=trim("SCALARS Y_SUM double")
      ! write(IOfile,'(a)') adjustl(trim(str))
      ! str=trim("LOOKUP_TABLE default")
      ! write(IOfile,'(a)') adjustl(trim(str))
      ! do i=1,NE
      !   write(str3,'(es30.12)') YSUM(i)
      !   write(IOfile,'(a)') adjustl(trim(str3))
      ! enddo

      do ks=1, Nvar-Nflow
        ! write(IOfile,*) 
        str=trim("SCALARS Production_rate_"//trim(getstring(ks))//" double")
        write(IOfile,'(a)') adjustl(trim(str))
        str=trim("LOOKUP_TABLE default")
        write(IOfile,'(a)') adjustl(trim(str))
        do i=1,NE
          write(str3,'(es30.12)') Production_rate(ks,i)
          write(IOfile,'(a)') adjustl(trim(str3))
        enddo
      enddo

      CLOSE(IOfile)

      if (myProc==1) then

        filename= trim("master"//trim(getstring(jout))//".visit")
        OPEN(IOfile,FILE=filename)

        str2= trim("!NBLOCKS "//trim(getstring(NPar)))
        write(IOfile,'(a)') adjustl(trim(str2))

        do i=1,NPar
          filename= trim("output"//trim(getstring(jout))//"_"//trim(getstring(i))//".vtk")
          write(IOfile,'(a)') adjustl(trim(filename)) 
        enddo

        close(IOfile)
      endif
      
    deallocate(RHO)
    deallocate(U)
    deallocate(V)
    deallocate(PRE)
    deallocate(MA)
    deallocate(TEM)
    deallocate(Y)
    ! deallocate(YSUM)
      
END subroutine

!!========================================================
!! output the start solution for ROM
!!--------------------------------------------------------
subroutine save_start_solution
    use mainvar
    use Parallel_Var
    implicit none
    integer,parameter:: IOfile=92
    integer i,k

    if(IfSteady) return
    if(IfRightDomain) return
    if(ttime < save_end_time - 0.1*dt .or. ttime > save_end_time + 0.1*dt) return
      
    open(unit=IOfile,file="sfvp"//trim(strProc)//".sav")
    write(IOfile,*) ttime
    write(IOfile,*) NE
    do i=1, NE
      write(IOfile,*) CellXY(1,i), CellXY(2,i), (PA(k,i),k=1,NVar)
    enddo
    close(IOfile)
    
end subroutine

!!======================================
!! back and recover data
!!--------------------------------------
subroutine BackupData
    use mainvar
    use Parallel_Var
    implicit none
    integer,parameter:: IOfile=92
    integer i,k

    ! if(.not. IfSteady) return
      
    open(unit=IOfile,file="fvp"//trim(strProc)//".sav")
    write(IOfile,*) ttime
    write(IOfile,*) NE
    do i=1, NE
      write(IOfile,*) CellXY(1,i), CellXY(2,i), (PA(k,i),k=1,NVar)
    enddo
    close(IOfile)
    
end subroutine
      
      
subroutine Read_Backup
    use mainvar
    use Parallel_Var
    implicit none
    double precision:: epsilon= 1.0e-12
    integer,parameter:: IOfile=92
    double precision:: temp(NVar),cx,cy
    integer i,j,k,ifile,nc
    character(len=50):: str
    integer,allocatable:: num_find_backup(:)
    character*50,external:: getstring

    allocate(num_find_backup(NE)); num_find_backup=0
      
    write(*,*) 'reading backup data...'

    if(IfRightDomain) then
      do ifile=1,NPreviousPar

        open(unit=IOfile,file="sfvp"//trim(getstring(ifile))//".sav", status="old",action="read")
        read(IOfile,*) ttime
        read(IOfile,*) nc
        do j=1, nc
          read(IOfile,*) cx,cy,(temp(k),k=1,NVar)
          do i=1,NE
            if(abs(CellXY(1,i)-cx)+abs(CellXY(2,i)-cy)<epsilon) then
              num_find_backup(i)=num_find_backup(i)+1
              PA(:,i)=temp(:)
            endif 
          enddo
        enddo
        close(IOfile)
      enddo

      if(any(1/=num_find_backup(:))) then
        write(*,*) 'Error in the finding the target cell in the backup file'
        call MPI_STOP
      endif
    else
      
      open(unit=IOfile,file="fvp"//trim(strProc)//".sav", status="old",action="read")
      read(IOfile,*) ttime
      read(IOfile,*) nc
      if(nc/=NE) then
        write(*,*) 'Error in the number of cells when reading backup'
        call MPI_STOP
      endif
      do i=1, NE
        read(IOfile,*) cx,cy,(PA(k,i),k=1,NVar)
      enddo
      close(IOfile)
    endif

    deallocate(num_find_backup)

    call Connect_CellAverage

    call primitve_to_conserved_variables

    write(*,*) 'finish reading backup solution ...'

    return
end subroutine
