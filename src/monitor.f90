!!========================================================
!! output the flow variable values of the monitor point
!!--------------------------------------------------------
subroutine write_monitor_point_value
    use mainvar
    implicit none
    integer,parameter:: IOfile= 115 
    double precision:: xy(2),dxyL(2),rof,uf,vf,pf,tf,vvh(MaxOrder)
    integer:: NL,id_point,ig
    integer,save :: isave=0
    character*50 :: filename
    character*50,external:: getstring

    if(num_monitor_points==0) return
    if(IfSteady) return

    do ig=1,num_monitor_points

      id_point= index_monitor_global(ig)

      if(IfRightDomain) then
        filename= trim("monitor_"//trim(getstring(id_point))//".dat")
      else
        filename= trim("monitor_fom_"//trim(getstring(id_point))//".dat")
      endif
    
      if(isave>0) then
        open(IOfile, file=filename, status="old", position="append", action="write")
      else
        open(IOfile, file=filename)
      endif

      !! compute the values of the monitor point
      ! xy(:)= (/x_monitor_point(id_point), y_monitor_point(id_point) /)
      NL= monitor_point_cell(id_point)
      
      ! dxyL(:)= xy(:)- CellXY(:,NL)
      ! call fbaseAll(dxyL,vvh,NL)

      ! write(IOfile,'(7e20.8)') ttime,mass_flow,rof,uf,vf,pf,tf
      ! write(IOfile,*) ttime, PA(:,NL)
      write(IOfile,'(5es20.8)') ttime,PA(1:4,NL)

      close(IOfile)
    enddo

    isave=isave + 1
    
    return
end subroutine

subroutine find_monitor_point
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer::i,num_find_points,num_find_points_sum,ierr,ig
    double precision:: xv1(2),xv2(2),xv3(2),xv4(2),xq(2),area1,area2,area3,area4
    double precision,external:: Triarea

    allocate(monitor_point_cell(num_global_monitor_points)); monitor_point_cell(:)= -1000
    allocate(index_monitor_global(num_global_monitor_points)); index_monitor_global= -1000

    num_monitor_points=0

    do ig=1,num_global_monitor_points

      xq=(/x_monitor_point(ig),y_monitor_point(ig)/)

      num_find_points=0

      do i=1,NE

        xv1(:)= Coor(:,n(1,i))
        xv2(:)= Coor(:,n(2,i))
        xv3(:)= Coor(:,n(3,i))
        xv4(:)= Coor(:,n(4,i))

        area1=Triarea(xv1,xv2,xq)
        area2=Triarea(xv2,xv3,xq)
        area3=Triarea(xv3,xv4,xq)
        area4=Triarea(xv4,xv1,xq)
    
        if(abs(area1+area2+area3+area4-vol(i))<1.e-12) then
          num_find_points= num_find_points+1
          monitor_point_cell(ig)= i
          num_monitor_points= num_monitor_points+1
          index_monitor_global(num_monitor_points)= ig
        endif
      enddo

      call MPI_AllReduce(num_find_points,num_find_points_sum, 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(num_find_points_sum>1) then
        write(*,*) 'find more than one monitor cells, num_find_points_sum=',num_find_points_sum, 'for point:',xq
        stop
      elseif(num_find_points_sum==1) then
        if(num_find_points==1) then
          write(*,*) 'monitor_point_cell=',monitor_point_cell(ig),'in processor',MyProc,'for point:',xq
          write(*,*) 'cell center of the monitor cell, x=', cellxy(1,monitor_point_cell(ig)),'y=',cellxy(2,monitor_point_cell(ig))
        endif
      else
        write(*,*) 'error: monitor cell for point',xq,'is not found'
        stop
      endif

    enddo
    
end subroutine

double precision function  Triarea(xv1,xv2,xv3)
     implicit none
     double precision :: xv1(2),xv2(2),xv3(2),x1,x2,x3,y1,y2,y3
   
    x1=xv1(1); y1= xv1(2)
    x2=xv2(1); y2= xv2(2)
    x3=xv3(1); y3= xv3(2)
    Triarea= 0.5*abs( (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1) )
     
end function