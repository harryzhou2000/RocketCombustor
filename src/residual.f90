!!==========================================================
!!                  Residual  
!!==========================================================
subroutine Residual(sumro,romax)
    use mainvar
    implicit none
    double precision::sumro(Nvar),romax
    integer:: i
    
    sumro=0.;  romax= 0.
!    sumvol=0.
    do i=1, NE
!      sumvol=sumvol+vol(i)
     sumro= sumro+ vol(i)*abs(hWA(:,i))
      ! sumro= sumro+ abs(hWA(1,i))
      romax= max(romax,WA(1,i))
    enddo
!    sumro=sumro/sumvol
    
    return
end subroutine

subroutine inner_iteration_residual(istep,inner_rk,RoRes)
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer, parameter:: iofile= 777
    integer,save:: isave=0
    double precision:: Res(Nvar),RoRes(1),Res_sum(Nvar),Romax(1),Romax_All(1)
    integer:: istep,inner_rk,ierr

    call Residual(Res, Romax(1))
    call MPI_AllReduce(Res,Res_sum,Nvar,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(Romax,Romax_All,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

    Res= Res_sum/total_volume
    Romax= Romax_All

    RoRes= sum(Res(Nflow+1:Nvar)) !Res(1)
    

    if( myProc==1 ) then

        filename= trim("converge_uns.dat")
        if(isave>0) then
            open(iofile, file=filename, status="old", position="append", action="write")
        else
            open(iofile, file=filename)
            WRITE(iofile,*)'VARIABLES = "Time step","Inner iteration step", "L1(Res)","Romax"'
        endif
        isave=isave + 1

        write(iofile,*) istep,inner_rk,Res,Romax
        close(iofile)
    endif
    
    return
end subroutine

subroutine step_residual(istep,cputime,RoRes)
    use mainvar
    use Parallel_Var
    implicit none
    include 'mpif.h'
    integer, parameter:: iofile= 666
    integer,save:: isave=0
    double precision:: Res(Nvar),Res_sum(Nvar),Romax(1),Romax_All(1),RoRes(1),cputime
    integer:: istep,ierr

    call Residual(Res,Romax(1))
    call MPI_AllReduce(Res,Res_sum,Nvar,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(Romax,Romax_All,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

    Res= Res_sum/total_volume
    Romax= Romax_All

    ! RoRes= Res(1)
    RoRes= sum(Res(Nflow+1:Nvar)) !Res(1)

    if( myProc==1 ) then

        filename= trim("converge.dat")
        if(isave>0) then
            open(iofile, file=filename, status="old", position="append", action="write")
        else
            open(iofile, file=filename)
            WRITE(iofile,*) 'VARIABLES = "Time step","CPU time (Sec)","L1(Res)","Romax"'
        endif
        isave=isave + 1

        write(iofile,*) istep,cputime,Res,Romax
        close(iofile)
    endif
    
    return
end subroutine