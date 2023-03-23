!!------------------------------------------------------
!! high order FVM code for 3D Navier-Stokes equs
!!   Wanai Li, 2011/11/18
!!----------------------------------
    
program main
    use mainvar
    implicit none
    
    !! geometry
    call input
    
    !! init reconstruction: form least-square matrix
    call Reconst_init
    
    call Partition
    
    !! output init flow-field
    !call NODE(0)
    
end program