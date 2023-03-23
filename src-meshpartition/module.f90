
module mainvar
    implicit none


    !! parallel info
    integer,save:: NPar
    integer,parameter :: MaxVir=10000,MaxCon=10

    integer,save:: MaxNei,Maxnei_bound

    integer,save :: Ninter,Nbound
    integer,allocatable:: itib(:),inter(:),bound(:),IsBoundCell(:)

    !! mesh size info
    integer,save:: NE,NNT,NF

     !! Mesh topological relation
    !! node
	integer,allocatable ::  NNodeNC(:),NodeNC(:,:)
    !!face
    integer,allocatable ::  FNode(:,:), NEIGHBOR(:,:), FProperty(:)
    integer,allocatable::   NodeFace(:,:),NumNFace(:)  ! ��¼����������棬��������
    !! cell
    integer,allocatable ::  N(:,:),NeighF(:,:), KFace(:),NeighC ( :,:), Neigh (:)

    !! Mesh geometry info
    real*8,parameter :: unit=1.  !!ע�⣬�����ٳ߶�unitӦ���Ǹ�������������
    !! node coordinates
    real*8,allocatable :: COOR(:,:)
    !! face
    real*8,allocatable :: VECX(:,:),SideOfs(:,:)
    !! cell 
    real*8,allocatable :: VOL(:),CELLXY(:,:),Ofs(:,:,:)
    
    !! curved boundary
    logical,parameter :: ExistCurveBoundary=.false.   !! ע�⣬������������ڵĳ������У�Ψһ��Ҫ�����Ĳ���
    integer,parameter :: WallID=-1
    !! curve and wall cell index
    integer,allocatable :: IsCurveCell(:),IfWallcell(:)
    !! curve face index
    integer,allocatable :: IsCurveFace(:)

    !! commen boundary information
    integer,allocatable:: IsCommenBoundaryNode(:,:)
    
!! connectivity
type Connect
   integer :: tarPar,tarCon, NListR, NListS
   integer :: ListR(0:2,MaxVir),ListS(2,MaxVir) !ListR(2,NListR),ListS(2,NListS)
   !! List(2,NList): 1: this partition local index(virtual); 2: target local index
end type Connect
    ! integer :: NConE(NPar)
    ! type(Connect),target :: ConE(MaxCon,NPar)  !! ConE(NConE,NPar)

    integer,allocatable :: NConE(:)
    type(Connect),target,allocatable :: ConE(:,:)  !! ConE(NConE,NPar)


type ConnectF
    integer :: tarPar,tarCon, NList, &
               List(MaxVir)
    !! List(:): face global index
end type ConnectF
    ! integer :: NConF(NPar)
    ! type(ConnectF),target :: ConF(MaxCon,NPar)

    integer,allocatable :: NConF(:)
    type(ConnectF),target,allocatable :: ConF(:,:)

end module
