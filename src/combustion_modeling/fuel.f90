!!=============================
!!  Fuel info
!!-----------------------------
module fuel
    implicit none

    double precision:: fuel_inlet_diameter, fuel_bottom, fuel_top
    double precision:: fuel_temp_in, fuel_pre_in, fuel_mass_flow
    double precision,allocatable:: fuel_mass_fraction(:)

end module

subroutine setup_fuel_info
    use species
    use fuel
    implicit none
    integer, parameter:: iofile=111
    character*50:: filename
    integer:: i

      call allocate_fuel_arrary

      filename= trim("combustion/fuel.in")
      OPEN(iofile,FILE=filename,STATUS='OLD')

      write(*,*) 'reading fuel info ...'

      read(iofile,*) 
      do i=1,Nsp
         read(iofile,*) fuel_mass_fraction(i)
      enddo

      read(iofile,*) 
      read(iofile,*) fuel_temp_in,fuel_pre_in

      read(iofile,*) 
      read(iofile,*) fuel_mass_flow

      close(iofile)
    
end subroutine

subroutine allocate_fuel_arrary
    use species
    use fuel
    implicit none

    allocate(fuel_mass_fraction(Nsp));fuel_mass_fraction=0.

    return
end subroutine


subroutine fuel_inlet_bc(pl,pr,ur,vr,tr,rr,Yr,Cpr,gamar,Rcpcvr)
    use species
    use fuel
    implicit none
    double precision:: pl,pr,ur,vr,tr,rr,Yr(Nsp),Cpr,gamar,Rcpcvr

    Yr(:)= fuel_mass_fraction(:)
    tr= fuel_temp_in
    call ComputeGasParameter(Yr,tr,Cpr,gamar,Rcpcvr)
    pr= pl
    rr= pr/(Rcpcvr*tr)
    ur= fuel_mass_flow/rr
    vr= 0.

    return
end subroutine
