!!=============================
!!  Oxidizer info
!!-----------------------------
module oxidizer
  implicit none

  double precision:: oxidizer_inlet_diameter, oxidizer_bottom, oxidizer_top
  double precision:: oxidizer_temp_in,oxidizer_pre_in,oxidizer_mass_flow
  double precision,allocatable:: oxidizer_mass_fraction(:)

end module

subroutine setup_oxidizer_info
    use species
    use oxidizer
    implicit none
    character*50:: filename
    integer, parameter:: iofile=123
    integer:: i

      call allocate_oxidizer_arrary

      filename= trim("combustion/oxidizer.in")
      OPEN(iofile,FILE=filename,STATUS='OLD')

      write(*,*) 'reading oxidizer info ...'

      read(iofile,*) 
      do i=1,Nsp
         read(iofile,*) oxidizer_mass_fraction(i)
      enddo

      read(iofile,*) 
      read(iofile,*) oxidizer_temp_in,oxidizer_pre_in

      read(iofile,*) 
      read(iofile,*) oxidizer_mass_flow

      close(iofile)
    
end subroutine

subroutine allocate_oxidizer_arrary
    use species
    use oxidizer
    implicit none

    allocate(oxidizer_mass_fraction(Nsp));oxidizer_mass_fraction=0.

    return
end subroutine

subroutine oxidizer_inlet_bc(pl,pr,ur,vr,tr,rr,Yr,Cpr,gamar,Rcpcvr)
    use species
    use oxidizer
    implicit none
    double precision:: pl,pr,ur,vr,tr,rr,Yr(Nsp),Cpr,gamar,Rcpcvr

    Yr(:)= oxidizer_mass_fraction(:)
    tr= oxidizer_temp_in
    call ComputeGasParameter(Yr,tr,Cpr,gamar,Rcpcvr)
    pr= pl
    rr= pr/(Rcpcvr*tr)
    ur= oxidizer_mass_flow/rr
    vr= 0.

    return
end subroutine

        