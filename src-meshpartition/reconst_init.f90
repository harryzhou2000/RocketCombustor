subroutine Reconst_init
    use mainvar
    integer :: NTemp(NE)
    
!! NeighF(3,NE)
    Ntemp(:)=0
    Nwall=0
    do i=1, NF
      NL = Neighbor(1,i)
      NR = Neighbor(2,i)
      
      if(NL<=0)stop
      Ntemp(NL)=Ntemp(NL)+1
      NeighF(Ntemp(NL),NL)= i
      if(NR>0 .and. Fproperty(i)/=100)then
        Ntemp(NR)=Ntemp(NR)+1
        NeighF(Ntemp(NR),NR)= i
      endif
    enddo
    
    
!! first neighbor cells: faces neighbored
    neigh(:) = 0
    Ninter=0
    NBound=0
    IsBoundCell=-100
    do i=1, NE
      do k=1,KFace(i)  !! 3
        iface = NeighF(k,i)    !! Î´ÓÐ¶¨Òå
        do kf=1,2
        ic= Neighbor(kf,iface)
        if( ic/=i .and. ic>0 )then
          Neigh(i)= Neigh(i)+1
          NeighC(Neigh(i),i)= ic
          
          Ofs(1:2,Neigh(i),i)= SideOfs(1:2,iface)
        endif
        enddo
      enddo
      
     if(Neigh(i)==KFace(i)) then
       Ninter=Ninter+1
       itib(i)=Ninter
       inter(Ninter)=i
     elseif(Neigh(i)>0 .and. Neigh(i)<KFace(i)) then
       IsBoundCell(i)=100
       Nbound=Nbound+1
       itib(i)=Nbound
       bound(nbound)=i
     else
        write(*,*) 'error: Neigh(i)> KFace(i) for cell', i
        stop
     endif
    enddo
    
    if(Nbound+Ninter/=NE) then
    write(*,*) 'error in boundcell number'
    stop
    endif 

    do i=1, NE
      NeighC (0,i)= i
      Ofs(1:2,0,i)= 0.
    enddo
    
    !!------------------------------------------
    write(*,*) 'reconstruction stencil,min,mincell/max,maxcell:',minval(Neigh),minloc(Neigh),&
                maxval(Neigh),maxloc(Neigh)
    
end subroutine


logical function contain(in,iarr,n)
    implicit none
    integer :: in,n,iarr(n),i
    do i=1,n
    if( in==iarr(i) )then
      contain= .true.
      return
    endif
    enddo
    contain=.false.
end function

