!!=============================
!! transport property
!!-----------------------------
module transport_property
    implicit none

    double precision:: transport_lowest_tem, transport_critical_tem, transport_highest_tem

    double precision,allocatable:: low_mu_A(:),low_mu_B(:),low_mu_C(:),low_mu_D(:),&
                                   high_mu_A(:),high_mu_B(:),high_mu_C(:),high_mu_D(:),&
                                   low_K_A(:),low_K_B(:),low_K_C(:),low_K_D(:),&
                                   high_K_A(:),high_K_B(:),high_K_C(:),high_K_D(:)

end module

subroutine setup_transport_property
    use species
    use transport_property
    implicit none
    integer:: i
    character*50:: filename
    integer,parameter :: iofile=115

    call allocate_transport_array

    filename= trim("combustion/transport.in")
    OPEN(iofile,FILE=filename,STATUS='OLD')

    write(*,*) 'reading transport property info ...'

    read(iofile,*)
    read(iofile,*) transport_lowest_tem,transport_critical_tem,transport_highest_tem

    ! robust control
    ! transport_lowest_tem = transport_lowest_tem - 1.

    read(iofile,*)
    do i=1,Nsp
        read(iofile,*) low_mu_A(i),low_mu_B(i),low_mu_C(i),low_mu_D(i)
        read(iofile,*) high_mu_A(i),high_mu_B(i),high_mu_C(i),high_mu_D(i)
        read(iofile,*) low_K_A(i),low_K_B(i),low_K_C(i),low_K_D(i)
        read(iofile,*) high_K_A(i),high_K_B(i),high_K_C(i),high_K_D(i)
    enddo

    close(iofile)

    return
end subroutine

subroutine allocate_transport_array
    use species
    use transport_property
    implicit none

    allocate(low_mu_A(Nsp));low_mu_A=0.
    allocate(low_mu_B(Nsp));low_mu_B=0.
    allocate(low_mu_C(Nsp));low_mu_C=0.
    allocate(low_mu_D(Nsp));low_mu_D=0.

    allocate(high_mu_A(Nsp));high_mu_A=0.
    allocate(high_mu_B(Nsp));high_mu_B=0.
    allocate(high_mu_C(Nsp));high_mu_C=0.
    allocate(high_mu_D(Nsp));high_mu_D=0.

    allocate(low_K_A(Nsp));low_K_A=0.
    allocate(low_K_B(Nsp));low_K_B=0.
    allocate(low_K_C(Nsp));low_K_C=0.
    allocate(low_K_D(Nsp));low_K_D=0.

    allocate(high_K_A(Nsp));high_K_A=0.
    allocate(high_K_B(Nsp));high_K_B=0.
    allocate(high_K_C(Nsp));high_K_C=0.
    allocate(high_K_D(Nsp));high_K_D=0.

    return
end subroutine

subroutine compute_transport_property(tem2,yy,vmul,akmu)
    use species
    use transport_property
    implicit none
    double precision:: tem2,tem,vmul,akmu,xx(Nsp),yy(Nsp),mul(Nsp),Kl(Nsp),phil(Nsp),ww
    integer:: i,j

    ! robust control
    tem= max(transport_lowest_tem,tem2)

    ww= 1./sum(yy(:)/mole_weight(:))
    xx(:)= ww*yy(:)/mole_weight(:)

    if(transport_lowest_tem <= tem .and. tem <= transport_critical_tem) then
        do i=1,Nsp
            mul(i)= exp(low_mu_A(i)*log(tem) + low_mu_B(i)/tem + low_mu_C(i)/tem**2 + low_mu_D(i))
            Kl(i)= exp(low_K_A(i)*log(tem) + low_K_B(i)/tem + low_K_C(i)/tem**2 + low_K_D(i))
        enddo
    elseif(transport_critical_tem < tem .and. tem <= transport_highest_tem) then
        do i=1,Nsp
            mul(i)= exp(high_mu_A(i)*log(tem) + high_mu_B(i)/tem + high_mu_C(i)/tem**2 + high_mu_D(i))
            Kl(i)= exp(high_K_A(i)*log(tem) + high_K_B(i)/tem + high_K_C(i)/tem**2 + high_K_D(i))
        enddo
    else
        write(*,*) 'error: temperature',tem,'is not in range:(',transport_lowest_tem,',',transport_highest_tem,')'
        ! call MPI_STOP
        stop
    endif

    !! units conversion
    mul= mul* 1.0e-7 !! from micropoise to Pa*s
    Kl = Kl * 1.0e-4 !! from microW/cm-K to W/m-K
    
    do i=1,Nsp

        phil(i)=0.
        do j=1,Nsp
            phil(i)= phil(i) + xx(j)*(1. + sqrt(mul(i)/mul(j))*(mole_weight(j)/mole_weight(i))**0.25)**2 &
                           /sqrt(1. + mole_weight(i)/mole_weight(j))
        enddo
        phil(i)= 0.25*sqrt(2.)*phil(i)
    enddo

    vmul= sum(xx(:)*mul(:)/phil(:))
    akmu= 0.5*(sum(xx(:)*Kl(:)) + 1./sum(xx(:)/Kl(:)))
    
    return
end subroutine
