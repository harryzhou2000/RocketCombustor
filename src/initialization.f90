subroutine initialization
    use mainvar
    use species
    use fuel
    use oxidizer
    implicit none
    integer:: i
    double precision, parameter:: length_step= 1.01e-02,diameter_chamber= 4.4e-02
    double precision, parameter:: reactant_temperature= 2000., ignition_temperature= 2000. !300. !2200.
    double precision,parameter:: reactant_mass_fraction(4)= (/ 0., 0.11750305997552021, 0.18849449204406366, &
        0.6940024479804162 /)
    ! double precision, parameter:: reactant_mass_fraction(4)= (/ 0., 0.42, 0., 0.58 /)
    double precision:: rr,uu,vv,tt,pp,Cp,gama,Rcpcv,YY(Nsp), xc, yc, oxidizer_velocity, fuel_velocity, &
                       reactant_mass_flow,reactant_velocity

    call compute_fuel_oxidizer_diameter

    reactant_mass_flow= fuel_mass_flow + oxidizer_mass_flow
    reactant_mass_flow= reactant_mass_flow/diameter_chamber

    oxidizer_mass_flow= oxidizer_mass_flow/oxidizer_inlet_diameter
    fuel_mass_flow= 0.5*fuel_mass_flow/fuel_inlet_diameter

    write(*,*) 'oxidzier mass flow= ', oxidizer_mass_flow
    write(*,*) 'fuel mass flow= ', fuel_mass_flow

    ttime=0.

    !! Attention here, we modify the fuel_mass_flow
    ! fuel_mass_flow= oxidizer_mass_flow !uu*rr
    p_static= p_out_ratio*oxidizer_pre_in  !!0.5,0.75
    p_out= p_static

    write(*,*) 'the static exit pressure is p_static=',p_static

    YY(:)= fuel_mass_fraction(:)
    tt= fuel_temp_in 
    call ComputeGasParameter(YY,tt,Cp,gama,Rcpcv)
    pp= p_out
    rr= pp/(Rcpcv*tt)
    fuel_velocity= fuel_mass_flow/rr

    write(*,*) 'fuel_velocity=',fuel_velocity


    YY(:)= oxidizer_mass_fraction(:)
    tt= oxidizer_temp_in 
    call ComputeGasParameter(YY,tt,Cp,gama,Rcpcv)
    pp= p_out
    rr= pp/(Rcpcv*tt)
    oxidizer_velocity= oxidizer_mass_flow/rr

    write(*,*) 'oxidizer_velocity=',oxidizer_velocity
       
    YY(:)= reactant_mass_fraction(:)
    tt= reactant_temperature
    call ComputeGasParameter(YY,tt,Cp,gama,Rcpcv)
    pp= p_out
    rr= pp/(Rcpcv*tt)
    reactant_velocity= reactant_mass_flow/rr

    write(*,*) 'reactant_velocity=',reactant_velocity

    !! special treatment
    ! fuel_velocity= reactant_velocity
    ! oxidizer_velocity= reactant_velocity

    do i=1,NE
        xc= cellxy(1,i)
        yc= cellxy(2,i)

        ! pp= oxidizer_pre_in + (p_out- oxidizer_pre_in)*(xc+0.04)/0.14
        if(xc<0) then
            if(yc> fuel_bottom ) then !! fuel injector
                YY(:)= fuel_mass_fraction(:)
                pp= p_out !fuel_pre_in
                tt= fuel_temp_in
                uu= fuel_velocity

                ! if(xc>= -length_step) tt= ignition_temperature 

                vv= 0.


            else !! oxidizer injector
               ! if(yc>oxidizer_top) then
                !    YY(:)= oxidizer_mass_fraction(:) + ( yc - oxidizer_top )/(fuel_bottom - oxidizer_top)&
                 !                                       *( fuel_mass_fraction(:) - oxidizer_mass_fraction(:) )
                !else
                    YY(:)= oxidizer_mass_fraction(:)
                !endif
                
                pp= p_out
                tt= oxidizer_temp_in

                if(yc>oxidizer_top) then
                    ! tt= reactant_temperature
                    ! uu= oxidizer_velocity + ( yc - oxidizer_top )/(fuel_bottom - oxidizer_top)&
                    !                                    *( fuel_velocity - oxidizer_velocity )
                    ! YY(:)= oxidizer_mass_fraction(:) + ( yc - oxidizer_top )/(fuel_bottom - oxidizer_top)&
                    !                                    *( fuel_mass_fraction(:) - oxidizer_mass_fraction(:) )

                    YY(:)= fuel_mass_fraction(:)
                    ! pp= p_out !fuel_pre_in
                    ! tt= fuel_temp_in !ignition_temperature !reactant_temperature !fuel_temp_in
                    ! tt= ignition_temperature + ( yc - oxidizer_top )/(fuel_bottom - oxidizer_top)&
                                                       ! *( fuel_temp_in - ignition_temperature )
                    tt= fuel_temp_in   !ignition_temperature
                    uu= fuel_velocity

                else
                    uu= oxidizer_velocity
                    if(xc>= -length_step) tt= ignition_temperature
                endif

                vv= 0.

            endif
        else   !! reactant 
            ! pp= oxidizer_pre_in
            YY(:)= reactant_mass_fraction(:)
            pp= p_out  
            tt= reactant_temperature
            uu= reactant_velocity
            vv= 0.
            
        endif
        PA(1,i)= pp
        PA(2,i)= uu
        PA(3,i)= vv
        PA(4,i)= tt
        PA(Nflow+1:Nvar,i)= YY(1:Nvar-Nflow)
    enddo

    call Connect_CellAverage

    call primitve_to_conserved_variables
    
    return
end subroutine


subroutine compute_fuel_oxidizer_diameter
    use fuel
    use oxidizer
    use mainvar
    implicit none
    include 'mpif.h'
    integer:: n1,n2,i,ierr
    double precision:: min_fuel_bc_y, max_fuel_bc_y, min_oxidizer_bc_y, max_oxidizer_bc_y
    double precision:: min_y(2),min_y_all(2),max_y(2),max_y_all(2)

    min_oxidizer_bc_y= 1.0e10; max_oxidizer_bc_y=  -1.0e10
    min_fuel_bc_y    = 1.0e10; max_fuel_bc_y    =  -1.0e10

    do i=1, NF
      n1= FNode(1,i); n2= FNode(2,i)
      if( Neighbor(2,i)==-3 )then !! oxidizer
        min_oxidizer_bc_y= min(min(coor(2,n1),coor(2,n2)),min_oxidizer_bc_y)
        max_oxidizer_bc_y= max(max(coor(2,n1),coor(2,n2)),max_oxidizer_bc_y)
      endif

      if( Neighbor(2,i)==-4 )then !! fuel
        min_fuel_bc_y= min(min(coor(2,n1),coor(2,n2)),min_fuel_bc_y)
        max_fuel_bc_y= max(max(coor(2,n1),coor(2,n2)),max_fuel_bc_y)
      endif
    enddo

    min_y=(/min_oxidizer_bc_y,min_fuel_bc_y/)
    max_y=(/max_oxidizer_bc_y,max_fuel_bc_y/)

    call MPI_AllReduce(min_y, min_y_all, 2, MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_AllReduce(max_y, max_y_all, 2, MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)

    min_oxidizer_bc_y= min_y_all(1); min_fuel_bc_y= min_y_all(2)
    max_oxidizer_bc_y= max_y_all(1); max_fuel_bc_y= max_y_all(2)

    fuel_bottom= min_fuel_bc_y ; fuel_top= max_fuel_bc_y
    oxidizer_bottom= min_oxidizer_bc_y ; oxidizer_top= max_oxidizer_bc_y

    oxidizer_inlet_diameter= 2*(max_oxidizer_bc_y - min_oxidizer_bc_y)
    fuel_inlet_diameter= max_fuel_bc_y - min_fuel_bc_y

    write(*,*) 'Oxidizer inlet diameter=', oxidizer_inlet_diameter
    write(*,*) 'Fuel inlet diameter=', fuel_inlet_diameter

    return
end subroutine
