subroutine allocate_ghostcell_arrary
    use mainvar
    implicit none

    allocate(Hostcell(Nghostcell)); Hostcell=0
    allocate(Hostface(Nghostcell)); Hostface=0
    allocate(Ghostcell(MaxNei,NE)); Ghostcell=0
    allocate(Num_ghostcell(NE)); Num_ghostcell=0
    allocate(GPA(Nvar,Nghostcell)); GPA=0.
    allocate(vol_Ghostcell(Nghostcell)); vol_Ghostcell=0.
    allocate(Coor_Ghostnode(2,4,Nghostcell)); Coor_Ghostnode=0. 
    allocate(CellXY_Ghostcell(2,Nghostcell)); CellXY_Ghostcell=0.
    allocate(BCtype_Ghostcell(Nghostcell)); BCtype_Ghostcell=0

    return
end subroutine


subroutine construct_ghostcells
    use mainvar
    implicit none
    real*8 :: sav1,sav2,sav,sav1n,sav2n,distance1,distance2
    real*8 :: xv1(2),xv2(2),xv3(2),xv4(2)
    integer:: i,j,k,iface,k1,k2,k3,k4,n1,n2,n3,n4,kf
    real*8,external :: GaussInteg_Quad

    Nghostcell=0
    do i=1,NE
        Nghostcell= Nghostcell + kface(i) - Neigh(i)
    enddo

    call allocate_ghostcell_arrary


    j= 0

	do i=1,NE

        Num_ghostcell(i)=0
    
    	if(Neigh(i)<kface(i)) then
    		do kf=1,kface(i)
                
                iface= neighF(kf,i)

        		if(Neighbor(2,iface)<0) then
                    
                    j= j + 1
                    Hostcell(j)= i
                    Hostface(j)= iface
                    BCtype_Ghostcell(j)= -Neighbor(2,iface)

                    Num_ghostcell(i)= Num_ghostcell(i) + 1
                    Ghostcell(Num_ghostcell(i),i)= j

                    n1= FNode(1,iface)
                    n2= FNode(2,iface)

                    do k1=1,kface(i)
                        if(n1==N(k1,i)) exit
                    enddo

                    k2= mod(k1,kface(i))+1

                    if(N(k2,i)/=n2) then
                        write(*,*) 'error in index setting of nodes when looking for ghostcells'
                        call MPI_STOP
                    endif

                    k3= mod(k2,kface(i))+1
                    k4= mod(k3,kface(i))+1

                    n3= N(k3,i)
                    n4= N(k4,i)


                    sav1= VECX(1,iface)
                    sav2= VECX(2,iface)
                    sav = max(sqrt(sav1**2+sav2**2),small)
                    sav1n=sav1/sav
                    sav2n=sav2/sav

                    Coor_Ghostnode(:,1,j)= Coor(:,n2)
                    Coor_Ghostnode(:,2,j)= Coor(:,n1)

                    distance1= (Coor(1,n1) - Coor(1,n4))*sav1n + (Coor(2,n1) - Coor(2,n4))*sav2n
                    if(distance1<0) then
                        write(*,*) 'error in calculating the distance when looking for ghost cells for node n1,n4'
                        write(*,*) 'distance=',distance1
                        write(*,*) 'n1,n2,n3,n4:',n1,n2,n3,n4
                        write(*,*) 'k1,k2,k3,k4:',k1,k2,k3,k4
                        write(*,*) 'x1,y1=',Coor(:,n1)
                        write(*,*) 'x2,y2=',Coor(:,n2)
                        write(*,*) 'x3,y3=',Coor(:,n3)
                        write(*,*) 'x4,y4=',Coor(:,n4)
                        write(*,*) 'BCtype=', BCtype_Ghostcell(j)
                        write(*,*) 'sav1,sav2=',sav1,sav2
                        write(*,*) 'sav1n,sav2n=',sav1n,sav2n
                        call MPI_STOP
                    endif
                    Coor_Ghostnode(:,3,j)= Coor(:,n4) + 2*distance1*(/sav1n,sav2n/)


                    distance2= (Coor(1,n2) - Coor(1,n3))*sav1n + (Coor(2,n2) - Coor(2,n3))*sav2n
                    if(distance2<0) then
                        write(*,*) 'error in calculating the distance when looking for ghost cells for node n2,n3'
                        write(*,*) 'distance=',distance2
                        write(*,*) 'n1,n2,n3,n4:',n1,n2,n3,n4
                        write(*,*) 'k1,k2,k3,k4:',k1,k2,k3,k4
                        write(*,*) 'x1,y1=',Coor(:,n1)
                        write(*,*) 'x2,y2=',Coor(:,n2)
                        write(*,*) 'x3,y3=',Coor(:,n3)
                        write(*,*) 'x4,y4=',Coor(:,n4)
                        write(*,*) 'BCtype=', BCtype_Ghostcell(j)
                        write(*,*) 'sav1,sav2=',sav1,sav2
                        write(*,*) 'sav1n,sav2n=',sav1n,sav2n
                        call MPI_STOP
                    endif
                    Coor_Ghostnode(:,4,j)= Coor(:,n3) + 2*distance2*(/sav1n,sav2n/)


                    xv1=Coor_Ghostnode(:,1,j)
                    xv2=Coor_Ghostnode(:,2,j)
                    xv3=Coor_Ghostnode(:,3,j)
                    xv4=Coor_Ghostnode(:,4,j)

                    vol_Ghostcell(j)= GaussInteg_Quad(fun, xv1,xv2,xv3,xv4)
                    CellXY_Ghostcell(:,j)= (/  GaussInteg_Quad(funx, xv1,xv2,xv3,xv4)/vol_Ghostcell(j), &
                                               GaussInteg_Quad(funy, xv1,xv2,xv3,xv4)/vol_Ghostcell(j) /)

                   
                    if( vol_Ghostcell(j)<0. )then
                    write(*,*) 'negative volume', j, vol_Ghostcell(j)
                    call MPI_STOP
                    endif      

                    !! print the information for ghost cell

                    ! write(*,*) 'n1,n2,n3,n4:',n1,n2,n3,n4
                    ! write(*,*) 'k1,k2,k3,k4:',k1,k2,k3,k4
                    ! write(*,*) 'x1,y1=',Coor(:,n1)
                    ! write(*,*) 'x2,y2=',Coor(:,n2)
                    ! write(*,*) 'x3,y3=',Coor(:,n3)
                    ! write(*,*) 'x4,y4=',Coor(:,n4)
                    ! write(*,*) 'BCtype=', BCtype_Ghostcell(j)
                    ! write(*,*) 'x1,y1=',xv1(:)
                    ! write(*,*) 'x2,y2=',xv2(:)
                    ! write(*,*) 'x3,y3=',xv3(:)
                    ! write(*,*) 'x4,y4=',xv4(:)

                    ! write(*,*) 'distance(n1,n4)=',distance1
                    ! write(*,*) 'distance(n2,n3)=',distance2
                    ! write(*,*) 'sav1n,sav2n=',sav1n,sav2n

                    ! write(*,*) 'cell center of host cell:',CellXY(:,i)
                    ! write(*,*) 'cell center of ghost cell:',CellXY_Ghostcell(:,j)
                    ! write(*,*) 'volume of ghost cell:',vol_Ghostcell(j)

            	endif
            enddo
    	endif

        if(Num_ghostcell(i)+Neigh(i)/=kface(i)) then
            write(*,*) 'error in number of ghost cells for cell:', i
            write(*,*) 'Neigh(i)=',Neigh(i),'Num_ghostcell(i)=',Num_ghostcell(i),'kface(i)=',kface(i)
            call MPI_STOP
        endif

    enddo

    if(j/=Nghostcell) then
        write(*,*) 'error in number of ghost cells'
        call MPI_STOP
    endif

contains
    real*8 function fun(xvec)
        implicit none
        real*8 :: xvec(2)
        fun= 1.
    end function

    real*8 function funx(xvec)
        implicit none
        real*8 :: xvec(2)
        funx= xvec(1)
    end function
    
    real*8 function funy(xvec)
        implicit none
        real*8 :: xvec(2)
        funy= xvec(2)
    end function

end subroutine
