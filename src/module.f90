
module mainvar
    implicit none

    logical,save:: If_turn_on_combustion
    
    !! number of flow variables, should be 4 + #turbulenceEquations
    integer,parameter:: NFlow=4

    !! total number of primitive/conserved variables
    integer, save:: NVar

    !! order setting
    integer, parameter:: MaxOrder=2    !! 2nd order
    integer, parameter:: NG=1, NG_f=2, NG_V=1, NG_V4=1

    integer,save:: MaxNei
    
    real*8,parameter :: small=1.0e-16, pi= 4*datan(1.d0)

    !! mesh size info
    integer,save:: NE,NE_Vir,NNT,NNT_vir,NF,NF_vir

    !! Mesh topological relation
    !! node
    integer,allocatable ::  NNodeNC(:),NodeNC(:,:)
    !!face
    integer,allocatable ::  FNode(:,:), NEIGHBOR(:,:), FProperty(:)
    integer,allocatable ::   NodeFace(:,:),NumNFace(:)  
    !! cell
    integer,allocatable ::  N(:,:),NeighF(:,:), KFace(:),NeighC ( :,:), Neigh (:)

    !! Mesh geometry info
    real*8, save :: unit
    !! total volume
    double precision:: total_volume
    !! node coordinates
    real*8,allocatable :: COOR(:,:)
    !! face
    real*8,allocatable :: VECX(:,:),SideOfs(:,:)
    !! cell 
    real*8,allocatable :: VOL(:),CELLXY(:,:), reflen(:,:),Ofs(:,:,:)

    !! Ghost cell array
    integer,save:: Nghostcell
    integer,allocatable:: Hostcell(:),Hostface(:),Ghostcell(:,:),Num_ghostcell(:),BCtype_Ghostcell(:)
    real*8,allocatable:: GPA(:,:),vol_Ghostcell(:), Coor_Ghostnode(:,:,:), &
                  CellXY_Ghostcell(:,:)

    !! geometry weight
    real*8,allocatable:: wdis(:,:)
    
    !! physical data, reconstruction coefficients,flux
    real*8,allocatable ::  WA(:,:), WP(:,:), WP2(:,:), PA(:,:),&
                           GradC(:,:,:), grad_tmp(:,:,:), &   
                           hWA(:,:), hWAP(:,:,:),Deltaeig(:), &
                           XYH(:,:), &
                           Smooth(:)

    double precision,allocatable:: production_rate(:,:), production_rate_sum(:), flux_difference(:), &
                                   flux_difference_nosource(:),sum_of_mass_fraction(:),phi_cell(:,:)
    
    !! spectral radius of cell interfaces
    real*8,allocatable :: alpaF(:)
    
    !! reconstruction matrices
    real*8,allocatable :: AMLS(:,:,:)
    
    !! local time step
    real*8,allocatable :: step(:)
    !! evolution parameter
    integer,save :: implicit_time_method
    real*8,save  :: cfl, ttime, dt, totaltime, save_start_time, save_end_time
    integer,save :: InLoop,InnerLoop,NScreenOutput,NFileOutput,NFileBackup,nstep_max,NWriteMonitor
    logical,save :: IfSteady,IfViscous,IfConstantViscousity,IfReadBackup,IfExplicit,IfFixDt
    real*8,save  :: relTol 
    
    !! for implicit RK 
    real*8,save:: asj(6,6),bs(6),cs(6)
    real*8,allocatable:: fixed(:,:)

    !! gaussian points and weigts    
    real*8,save :: wxx(NG),wgg(NG),wxxV(3,NG_V),wggV(NG_V),wxxV4(2,NG_V4),wggV4(NG_V4)
    
    !! wall boundary
    integer,parameter :: WallID=-1

    !! limiting info
    integer, save:: limit_ID, n_wbap
    real*8, save:: rc_detector

    !! common boundary information
    integer,save:: num_common_boundary
    integer,allocatable:: IsCommonBoundaryNode(:,:),IsCommonBoundaryFace(:)

    !! Domain splitting information
    logical,save:: IfRightDomain
    integer,save:: NPreviousPar
    integer,allocatable:: index_common_boundary(:,:)

    !! monitor point
    integer,save:: num_global_monitor_points,num_monitor_points
    integer,allocatable:: monitor_point_cell(:),index_monitor_global(:)
    double precision,allocatable:: x_monitor_point(:),y_monitor_point(:)

    double precision,save:: p_out,p_out_ratio,p_static
    double precision,save:: A_p, f_p

contains
subroutine fbaseAll(xyz,vvh,icell)
    implicit none
    integer :: icell
    real*8 :: xyz(2),vvh(MaxOrder), x,y   !! xy(:): must be local coordindate
    
    x= xyz(1)/reflen(1,icell)
    y= xyz(2)/reflen(2,icell)
    
!    !! 2nd-order
    vvh(:)=(/ x,y /)
    
    return
end subroutine

subroutine fbaseGrad_dxdy(xyz,vvh,icell)
    implicit none
    integer :: icell
    real*8 :: xyz(2),vvh(MaxOrder,2), dx,dy   !! xy(:): must be local coordindate
   
    dx=1./reflen(1,icell)
    dy=1./reflen(2,icell)
    
!    !! 4th order
    vvh(:,1)=(/ 1.,0./)*dx
    vvh(:,2)=(/ 0.,1./)*dy

    return
end subroutine

function limited_mass_fraction(x)
    implicit none
    double precision, intent(in):: x(:)
    double precision:: limited_mass_fraction(size(x))
    double precision:: sum_x
    
    !limited_mass_fraction= x

    limited_mass_fraction= max(0.,x)
    limited_mass_fraction= min(1.,limited_mass_fraction)

    sum_x= sum(limited_mass_fraction)

    limited_mass_fraction= limited_mass_fraction/sum_x
    
end function

end module

function getstring(jout) result(str)
    implicit none
    integer:: jout
    character*50:: str

    if(jout<10)then
        write(str,"(i1)") jout
      elseif(jout>=10 .and. jout<100) then
        write(str,"(i2)") jout
      elseif(jout>=100 .and. jout<1000) then
        write(str,"(i3)") jout
      elseif(jout>=1000 .and. jout<10000) then
        write(str,"(i4)") jout
      elseif(jout>=10000 .and. jout<100000) then
        write(str,"(i5)") jout
      elseif(jout>=100000 .and. jout<1000000) then
        write(str,"(i6)") jout
      elseif(jout>=1000000 .and. jout<10000000) then
        write(str,"(i7)") jout
      elseif(jout>=10000000 .and. jout<100000000) then
        write(str,"(i8)") jout
      else
        write(*,*) 'Error: time step index', jout, 'larger than 10**8 when writting .visit file'
        call MPI_STOP
      endif
    
end function
