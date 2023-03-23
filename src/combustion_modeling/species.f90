!!=========================
!! Species info
!!-------------------------
module species
    implicit none

    integer:: Nsp
    double precision:: Rconstant
    double precision,allocatable:: mole_weight(:),R_species(:)
    character*50,allocatable:: species_name(:)

end module

subroutine setup_species_info
    use species
    implicit none
    integer, parameter:: iofile=113
    integer:: i
    character*50:: filename


      filename= trim("combustion/species.in")
      OPEN(iofile,FILE=filename,STATUS='OLD')

      write(*,*) 'reading species info ...'

      read(iofile,*) 
      read(iofile,*) Nsp

      call allocate_species_arrary

      read(iofile,*)
      do i=1,Nsp
        read(iofile,*) species_name(i),mole_weight(i)
      enddo

      read(iofile,*)
      read(iofile,*) Rconstant

      close(iofile)

      do i=1,Nsp
        R_species(i)= Rconstant/mole_weight(i)
      enddo

      write(*,*) 'Species:'
      do i=1,Nsp
        write(*,*) trim(species_name(i)), mole_weight(i), R_species(i)
      enddo
    
end subroutine

subroutine setup_number_of_variables
    use mainvar
    use species
    implicit none

    Nvar= Nflow + Nsp - 1

    write(*,*) 'Important reminder:'
    write(*,*) 'We set the #MassFractions = #Species - 1'
    write(*,*) 'Because the last species is a diluted gas'

    return
end subroutine

subroutine allocate_species_arrary
    use species
    implicit none

    allocate(species_name(Nsp));species_name=''

    allocate(mole_weight(Nsp));mole_weight=0.

    allocate(R_species(Nsp));R_species=0.

    return
end subroutine