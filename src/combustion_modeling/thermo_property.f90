!!=============================
!! thermo property
!!-----------------------------
module thermo_property
    implicit none

    double precision:: Tref

    integer:: order_of_cp_polynomial

    double precision:: thermo_lowest_tem,thermo_critical_tem,thermo_highest_tem

    double precision,allocatable:: hof_species(:),low_coef_Cp_species(:,:),&
                                   high_coef_Cp_species(:,:), low_enthalpy_b1(:), high_enthalpy_b1(:)

end module

subroutine setup_thermo_property
    use species
    use thermo_property
    implicit none
    character*50:: filename
    integer, parameter:: iofile=221
    integer:: i

    filename= trim("combustion/thermo.in")
    OPEN(iofile,FILE=filename,STATUS='OLD')

    write(*,*) 'reading thermo property info ...'

    read(iofile,*)
    read(iofile,*) Tref

    read(iofile,*) 
    read(iofile,*) order_of_cp_polynomial

    call allocate_thermo_arrary

    read(iofile,*) 
    read(iofile,*) thermo_lowest_tem,thermo_critical_tem,thermo_highest_tem
    do i=1,Nsp
        read(iofile,*) low_coef_Cp_species(1:order_of_cp_polynomial,i)
        read(iofile,*) low_enthalpy_b1(i)
        read(iofile,*) high_coef_Cp_species(1:order_of_cp_polynomial,i)
        read(iofile,*) high_enthalpy_b1(i)
    enddo

    read(iofile,*) 
    do i=1,Nsp
        read(iofile,*) hof_species(i)
        hof_species(i)= hof_species(i)/mole_weight(i)
        write(*,*) 'Heat of Formation:'
        write(*,*) trim(species_name(i)), hof_species(i)
    enddo

    close(iofile)

    return
end subroutine

subroutine allocate_thermo_arrary
    use species
    use thermo_property
    implicit none

    allocate(hof_species(Nsp));hof_species=0.

    allocate(low_coef_Cp_species(order_of_cp_polynomial,Nsp));low_coef_Cp_species=0.
    allocate(high_coef_Cp_species(order_of_cp_polynomial,Nsp));high_coef_Cp_species=0.

    allocate(low_enthalpy_b1(Nsp));low_enthalpy_b1=0.
    allocate(high_enthalpy_b1(Nsp));high_enthalpy_b1=0.

    return
end subroutine


subroutine ComputeGasParameter(Y, tem, Cp, gama, Rcpcv)
    use species
    use thermo_property
    implicit none
    integer:: i,j
    double precision :: Y(Nsp), tem, Cp_species(Nsp), Cp, gama, Rcpcv, cp_basis(order_of_cp_polynomial)

    do j=1,order_of_cp_polynomial
        cp_basis(j)=tem**(j-3)
    enddo

    if(thermo_lowest_tem <= tem .and. tem <= thermo_critical_tem) then
        do i=1,Nsp
            Cp_species(i)= sum(low_coef_Cp_species(:,i)*cp_basis(:))*R_species(i)
        enddo
    elseif(thermo_critical_tem < tem .and. tem <= thermo_highest_tem) then
        do i=1,Nsp
            Cp_species(i)= sum(high_coef_Cp_species(:,i)*cp_basis(:))*R_species(i)
        enddo
    else
        write(*,*) 'error in ComputeGasParameter: temperature',tem,'is not in range:(',thermo_lowest_tem,',',thermo_highest_tem,')'
        ! call MPI_STOP
        stop
    endif

    Cp   = sum(Cp_species(:)*Y(:))

    ! write(*,*) "CP" , Cp

    Rcpcv= sum(R_species(:)*Y(:))

    gama = Cp/(Cp-Rcpcv) 
    
    return
end subroutine

subroutine compute_enthalpy(tem,entha)
    use species
    use thermo_property
    implicit none
    integer:: i,j
    double precision:: tem,entha(Nsp),enthalpy_basis(order_of_cp_polynomial),ref_basis(order_of_cp_polynomial)
 

    enthalpy_basis(1:2)=(/ -tem**(-2),log(tem)/tem/)
    do j=3,order_of_cp_polynomial
        enthalpy_basis(j)= tem**(j-3)/(j-2)
    enddo

    !! enthalpy at reference temperature equals to the heat of formation

    if(thermo_lowest_tem <= tem .and. tem <= thermo_critical_tem) then
        do i=1,Nsp
            entha(i)= R_species(i)*( tem*sum(low_coef_Cp_species(:,i)*enthalpy_basis(:)) + low_enthalpy_b1(i) )
        enddo
    elseif(thermo_critical_tem < tem .and. tem <= thermo_highest_tem) then
        do i=1,Nsp
            entha(i)= R_species(i)*( tem*sum(high_coef_Cp_species(:,i)*enthalpy_basis(:)) + high_enthalpy_b1(i) )
        enddo
    else
        write(*,*) 'error in compute_enthalpy: temperature',tem,'is not in range:(',thermo_lowest_tem,',',thermo_highest_tem,')'
        ! call MPI_STOP
        stop
    endif

    return
end subroutine

subroutine Newton_temperature(et,Y,t,t0,Cp,gama,Rcpcv)
    use species
    use thermo_property
    implicit none
    integer, parameter:: max_iter=200
    integer:: num_iter
    double precision:: et,Y(Nsp),t, t0, h, ft, dft, dtt, entha(Nsp), Cp, gama, Rcpcv
    double precision, parameter:: epsilon= 1.0e-5

    num_iter=0; t= t0 !!initial guess

    do while (num_iter<max_iter)

        num_iter= num_iter + 1
        
        ! Newton procedure
        call compute_enthalpy(t,entha)
        call ComputeGasParameter(Y, t, Cp, gama, Rcpcv)
        h  = sum(Y(:)*entha(:))
        ft = h  - Rcpcv*t - et
        dft= Cp - Rcpcv
        dtt= - ft/dft
        t= t + dtt
        
        t=max(t,thermo_lowest_tem)
        t=min(t,thermo_highest_tem)

        ! write(*,*) num_iter, t, ft, dft, dtt
        

        if(abs(dtt)<epsilon) exit
    enddo

    return
end subroutine


