
module neuralnetworks
    implicit none

    integer,save :: Ncommon_bound

    type lstm
        integer:: dim_x, dim_hidden_state
        double precision, allocatable:: weight(:,:), recurrent_weight(:,:), bias(:)
        double precision, allocatable:: x(:), hidden_state(:), cell_state(:)
    end type lstm

    type rnn
        integer:: n_in, n_out, n_steps
        double precision, allocatable:: weight_output(:,:), bias_output(:)
        double precision, allocatable:: mean_input(:), mean_output(:), deviation_input(:), deviation_output(:)
        double precision, allocatable:: input(:,:), output(:)

        type(lstm):: lstm_unit
    end type rnn

    !! boundary info
    type common_bound
        integer:: n_cells, n_dof, n_basis, gap_saved_solutions
        double precision,allocatable:: var_ref(:)
        double precision,allocatable:: boundary_value(:), basis(:,:)
        type(rnn):: network

        double precision, allocatable:: start_solution(:,:), right_solution(:), reduced_right_solution(:)

        !! solution points on the current processor
        integer:: num_cells_on_partition
        integer,allocatable:: order_of_solution_components(:), total_order_of_solution_components(:)
        double precision,allocatable:: soluton_on_partition(:)


        !!parallel gathering info
        integer, allocatable::RECVCOUNTS(:), DISPLS(:)

    end type common_bound

    type(common_bound),allocatable,target :: Common_bound_targets(:)


contains

function relu(x)
    implicit none
    integer:: n,i
    double precision, intent(in):: x(:)
    double precision:: relu(size(x))

    relu= max(0.,x)

end function

function swish(x)
    implicit none
    double precision,intent(in):: x(:)
    double precision:: swish(size(x))

    swish= x/(1+exp(-x))

end function

function hard_sigmoid(x)
    implicit none
    integer:: i
    double precision,intent(in):: x(:)
    double precision:: hard_sigmoid(size(x))

    do i=1,size(x)

        if(x(i) < -2.5) then
            hard_sigmoid(i)= 0.
        elseif(x(i) > 2.5) then
            hard_sigmoid(i)= 1.
        else
            hard_sigmoid(i)= 0.2*x(i) + 0.5
        endif
        
    enddo

end function

function softmax(x)
    implicit none
    double precision,intent(in):: x(:)
    double precision:: softmax(size(x)),exp_sum

    exp_sum= sum(exp(x))
    
    softmax= exp(x)/exp_sum
    
end function

function scaling_input(x,mean_input,deviation_input)
    implicit none
    double precision,intent(in):: x(:),mean_input(:),deviation_input(:)
    double precision:: scaling_input(size(x))

    scaling_input= (x-mean_input)/deviation_input
    
end function

function scaling_output(y,mean_output,deviation_output)
    implicit none
    double precision,intent(in):: y(:),mean_output(:),deviation_output(:)
    double precision:: scaling_output(size(y))

   
    scaling_output= y*deviation_output + mean_output
    
end function

subroutine evaluate_lstm(lstm_unit)
    implicit none
    type(lstm):: lstm_unit
    double precision:: vec(4*lstm_unit%dim_hidden_state)
    integer::ndim
    double precision, dimension(lstm_unit%dim_hidden_state):: it, ft, gt, ot

    vec= matmul(lstm_unit%x,lstm_unit%weight) + matmul(lstm_unit%hidden_state, lstm_unit%recurrent_weight) + lstm_unit%bias

    ndim= lstm_unit%dim_hidden_state
    it= hard_sigmoid(vec(1:ndim))
    ft= hard_sigmoid(vec(ndim+1:2*ndim))
    gt= tanh(vec(2*ndim+1:3*ndim))
    ot= hard_sigmoid(vec(3*ndim+1:4*ndim))

    !! cell & hidden states
    lstm_unit%cell_state= ft*lstm_unit%cell_state + it*gt
    lstm_unit%hidden_state= ot*tanh(lstm_unit%cell_state)

end subroutine

subroutine forward_prop(network)
    implicit none 
    type(rnn):: network
    integer:: i

    network%lstm_unit%hidden_state= 0.; network%lstm_unit%cell_state=0.

    do i=1,network%n_steps
        network%lstm_unit%x= scaling_input(network%input(:,i),network%mean_input,network%deviation_input)
        call evaluate_lstm(network%lstm_unit)
    enddo

    network%output=matmul(network%lstm_unit%hidden_state,network%weight_output) + network%bias_output

    network%output= scaling_output(network%output, network%mean_output, network%deviation_output)

end subroutine
    
end module

subroutine load_neural_network
    use mainvar
    use neuralnetworks
    implicit none
    integer,parameter:: iofile=10
    character*50 :: filename,str
    integer:: i,j,ia,ib
    type(common_bound),pointer:: tcommon_bound
    type(rnn),pointer:: tnetwork
    character*50,external:: getstring

    write(*,*) 'reading common boundary information'

    filename=trim("../COMMON_BOUND/info_common_boundary.dat")
    OPEN(iofile,FILE=filename,STATUS='OLD')
    read(iofile,*) Ncommon_bound
    allocate(Common_bound_targets(Ncommon_bound))
    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)
        read(iofile,*) j, tcommon_bound%n_cells, tcommon_bound%n_basis
        tcommon_bound%n_dof= tcommon_bound%n_cells*Nvar
        allocate(tcommon_bound%var_ref(tcommon_bound%n_dof));tcommon_bound%var_ref=0.
        allocate(tcommon_bound%boundary_value(tcommon_bound%n_dof));tcommon_bound%boundary_value=0.
        allocate(tcommon_bound%basis(tcommon_bound%n_dof,tcommon_bound%n_basis));tcommon_bound%basis=0.
    enddo
    close(iofile)
    
    write(*,*) 'reading basis, weights & bias ...'
    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)

        !! read reference solution
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/ref_solution_"//trim(getstring(i))//".dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) (tcommon_bound%var_ref(ib),ib=1,tcommon_bound%n_dof)
        close(iofile)

        !! read the reduced basis functions
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/basis_"//trim(getstring(i))//".dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) ((tcommon_bound%basis(ia,ib),ib=1,tcommon_bound%n_basis),ia=1,tcommon_bound%n_dof)
        close(iofile)


        !! Now read network stuff
        tnetwork => tcommon_bound%network

        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/info.dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) tnetwork%n_steps, tnetwork%n_in, tnetwork%n_out, tnetwork%lstm_unit%dim_hidden_state
        close(iofile)

        !! input and output of the network
        allocate(tnetwork%input(tnetwork%n_in,tnetwork%n_steps));tnetwork%input=0.
        allocate(tnetwork%output(tnetwork%n_out));tnetwork%output=0.

        !! mean & variance of the input & output
        allocate(tnetwork%mean_input(tnetwork%n_in));tnetwork%mean_input=0.
        allocate(tnetwork%mean_output(tnetwork%n_out));tnetwork%mean_output=0.
        allocate(tnetwork%deviation_input(tnetwork%n_in));tnetwork%deviation_input=0.
        allocate(tnetwork%deviation_output(tnetwork%n_out));tnetwork%deviation_output=0.

        !! mean & variance of input
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/mean_input.dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) (tnetwork%mean_input(ib),ib=1,tnetwork%n_in)
        close(iofile)

        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/deviation_input.dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) (tnetwork%deviation_input(ib),ib=1,tnetwork%n_in)
        close(iofile)

        !! mean & variance of output
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/mean_output.dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) (tnetwork%mean_output(ib),ib=1,tnetwork%n_out)
        close(iofile)

        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/deviation_output.dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) (tnetwork%deviation_output(ib),ib=1,tnetwork%n_out)
        close(iofile)

        !! allocate and load weights & bias of the lstm unit
        tnetwork%lstm_unit%dim_x= tnetwork%n_in

        allocate(tnetwork%lstm_unit%x(tnetwork%lstm_unit%dim_x)); tnetwork%lstm_unit%x=0.
        allocate(tnetwork%lstm_unit%hidden_state(tnetwork%lstm_unit%dim_hidden_state)); tnetwork%lstm_unit%hidden_state=0.
        allocate(tnetwork%lstm_unit%cell_state(tnetwork%lstm_unit%dim_hidden_state)); tnetwork%lstm_unit%cell_state=0.

        allocate(tnetwork%lstm_unit%weight(tnetwork%lstm_unit%dim_x,4*tnetwork%lstm_unit%dim_hidden_state))
        tnetwork%lstm_unit%weight=0.

        allocate(tnetwork%lstm_unit%recurrent_weight(tnetwork%lstm_unit%dim_hidden_state,4*tnetwork%lstm_unit%dim_hidden_state))
        tnetwork%lstm_unit%recurrent_weight=0.

        allocate(tnetwork%lstm_unit%bias(4*tnetwork%lstm_unit%dim_hidden_state))
        tnetwork%lstm_unit%bias=0.

        !! read weight
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/W.txt")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) ((tnetwork%lstm_unit%weight(ia,ib),ib=1,4*tnetwork%lstm_unit%dim_hidden_state),ia=1,tnetwork%lstm_unit%dim_x)
        close(iofile)

        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/R.txt")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) ((tnetwork%lstm_unit%recurrent_weight(ia,ib),ib=1,4*tnetwork%lstm_unit%dim_hidden_state),&
            ia=1,tnetwork%lstm_unit%dim_hidden_state)
        close(iofile)

        !! read bias
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/b.txt")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) (tnetwork%lstm_unit%bias(ib),ib=1,4*tnetwork%lstm_unit%dim_hidden_state)
        close(iofile)


        !! allocate and read output layer weight & bias
        allocate(tnetwork%weight_output(tnetwork%lstm_unit%dim_hidden_state,tnetwork%n_out)); tnetwork%weight_output=0.
        allocate(tnetwork%bias_output(tnetwork%n_out)); tnetwork%bias_output=0.

        !! read weight
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/Wo.txt")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) ((tnetwork%weight_output(ia,ib),ib=1,tnetwork%n_out),ia=1,tnetwork%lstm_unit%dim_hidden_state)
        close(iofile)

        !! read bias
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/bo.txt")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) (tnetwork%bias_output(ib),ib=1,tnetwork%n_out)
        close(iofile)

        !! allocate and read start-up reduced solutions
        filename=trim("../COMMON_BOUND/NET_"//trim(getstring(i))//"/start_solution_"//trim(getstring(i))//".dat")
        open(iofile,file=filename,STATUS='old')
        read(iofile,*) tcommon_bound%gap_saved_solutions
        allocate(tcommon_bound%start_solution(tnetwork%n_in,(tnetwork%n_steps-1)*tcommon_bound%gap_saved_solutions+1))
        tcommon_bound%start_solution=0.
        read(iofile,*) ((tcommon_bound%start_solution(ia,ib),ia=1,tnetwork%n_in),&
            ib=1,(tnetwork%n_steps-1)*tcommon_bound%gap_saved_solutions+1)
        close(iofile)

        !! right solution
        allocate(tcommon_bound%right_solution(tcommon_bound%n_dof));tcommon_bound%right_solution=0.
        allocate(tcommon_bound%reduced_right_solution(tcommon_bound%n_basis));tcommon_bound%reduced_right_solution=0.
    enddo

    call set_MPI_gatherv

    return
end subroutine

subroutine set_MPI_gatherv
    use mainvar
    use neuralnetworks
    use Parallel_Var
    implicit none
    include 'mpif.h'
    type(common_bound),pointer:: tcommon_bound
    integer:: i,j,NL,NR,n1,n2,id_face,id_face_number,icount_face,ierr,send_count


    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)
        tcommon_bound%num_cells_on_partition=0
    enddo


    do i=1,NF
        NR= Neighbor(2,i)
        if(NR<-1000) then
            id_face= index_common_boundary(1,i)
            tcommon_bound => Common_bound_targets(id_face)
            tcommon_bound%num_cells_on_partition= tcommon_bound%num_cells_on_partition + 1
        endif
    enddo

    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)
        allocate(tcommon_bound%order_of_solution_components(Nvar*tcommon_bound%num_cells_on_partition))
        tcommon_bound%order_of_solution_components=0

        allocate(tcommon_bound%soluton_on_partition(Nvar*tcommon_bound%num_cells_on_partition))
        tcommon_bound%soluton_on_partition=0.

        allocate(tcommon_bound%total_order_of_solution_components(tcommon_bound%n_dof))
        tcommon_bound%total_order_of_solution_components=0

        !! parallel info
        allocate(tcommon_bound%RECVCOUNTS(Npar)); tcommon_bound%RECVCOUNTS=0
        allocate(tcommon_bound%DISPLS(Npar)); tcommon_bound%DISPLS=0
    enddo

    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)
        tcommon_bound%num_cells_on_partition=0
    enddo

    do i=1,NF
        NR= Neighbor(2,i)

        if(NR<-1000) then

            id_face= index_common_boundary(1,i)
            id_face_number= index_common_boundary(2,i)

            tcommon_bound => Common_bound_targets(id_face)
            tcommon_bound%num_cells_on_partition= tcommon_bound%num_cells_on_partition + 1

            icount_face= tcommon_bound%num_cells_on_partition

            do j=1,Nvar
                tcommon_bound%order_of_solution_components(Nvar*(icount_face-1)+j)= Nvar*(id_face_number-1)+j
            enddo

        endif
    enddo



    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)

        send_count= Nvar*tcommon_bound%num_cells_on_partition

        call MPI_ALLGATHER(send_count,1,MPI_INT,tcommon_bound%RECVCOUNTS,1,MPI_INT,MPI_COMM_WORLD,ierr)

        

        do j=2,Npar
            tcommon_bound%DISPLS(j)= sum(tcommon_bound%RECVCOUNTS(1:j-1))
        enddo

        ! write(*,*) 'processor', MyProc, 'num_cells_on_partition', tcommon_bound%num_cells_on_partition
        ! write(*,*) 'processor', MyProc, 'order_of_solution_components', tcommon_bound%order_of_solution_components

        

        call MPI_ALLGATHERV(tcommon_bound%order_of_solution_components, &
            Nvar*tcommon_bound%num_cells_on_partition, &
            MPI_INT, tcommon_bound%total_order_of_solution_components, tcommon_bound%RECVCOUNTS,&
            tcommon_bound%DISPLS, MPI_INT,MPI_COMM_WORLD,ierr)

        ! if(MyProc==1) then
        !     write(*,*) tcommon_bound%RECVCOUNTS
        !     write(*,*) tcommon_bound%DISPLS
        !     write(*,*) tcommon_bound%total_order_of_solution_components
        ! endif
    enddo

    ! call MPI_STOP

    return
end subroutine

subroutine update_right_solution
    use mainvar
    use neuralnetworks
    use Parallel_Var
    implicit none
    include 'mpif.h'
    type(common_bound),pointer:: tcommon_bound
    integer:: i,j,ndim,NL,NR,n1,n2,id_face,id_face_number,icount_face,ierr
    double precision:: dxyL(2),xy(2),vvh(MaxOrder),sol(Nvar)
    double precision,allocatable:: RECVBUF(:)

    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)
        tcommon_bound%num_cells_on_partition=0
    enddo

    ! icount_face=0
    do i=1,NF
        NL= Neighbor(1,i)
        NR= Neighbor(2,i)

        if(NR<-1000) then
            n1= FNode(1,i)
            n2= FNode(2,i)
            xy(:)= 0.5*(Coor(:,n1)+Coor(:,n2))

            id_face= index_common_boundary(1,i)
            id_face_number= index_common_boundary(2,i)

            tcommon_bound => Common_bound_targets(id_face)
            tcommon_bound%num_cells_on_partition= tcommon_bound%num_cells_on_partition + 1
            icount_face= tcommon_bound%num_cells_on_partition

            !! Left and Right state
            dxyL(:)= xy(:)- CellXY(:,NL)
            call fbaseAll(dxyL,vvh,NL)

            do j=1,Nvar
                sol(j)= PA(j,NL) + sum(gradC(:,j,NL)*vvh(:))
            enddo

            tcommon_bound%soluton_on_partition(Nvar*(icount_face-1)+1:Nvar*icount_face)= sol(:)

        endif
    enddo

    ! write(*,*) 'there are', icount_face, 'boundary faces in processor', MyProc

    !! communication to get the entire right solution
    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)

        allocate(RECVBUF(tcommon_bound%n_dof))
        ! write(*,*) 'processor', MyProc, 'num_cells_on_partition', tcommon_bound%num_cells_on_partition
        ! write(*,*) 'processor', MyProc, 'order_of_solution_components', tcommon_bound%order_of_solution_components
        ! write(*,*) 'processor', MyProc, 'soluton_on_partition', tcommon_bound%soluton_on_partition

        ! call MPI_ALLGATHERV(tcommon_bound%soluton_on_partition, &
        !     Nvar*tcommon_bound%num_cells_on_partition, &
        !         MPI_DOUBLE_PRECISION, RECVBUF, tcommon_bound%RECVCOUNTS,&
        !         tcommon_bound%DISPLS, MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)

        write(*,*) 'Here we need to uncomment the MPI_ALLGATHERV line'
        call MPI_STOP
        

        do j=1,tcommon_bound%n_dof
            tcommon_bound%right_solution(tcommon_bound%total_order_of_solution_components(j))= RECVBUF(j)
        enddo

        ! if(MyProc==1) then
        ! write(*,*) 'processor', MyProc, 'RECVBUF', RECVBUF
        ! write(*,*) 'processor', MyProc, 'right_solution', tcommon_bound%right_solution
        ! endif

        deallocate(RECVBUF)
    enddo


    !! project the right solution onto reduced space
    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)

        !! compute reduced right solution
        tcommon_bound%reduced_right_solution= matmul(transpose(tcommon_bound%basis),&
            tcommon_bound%right_solution - tcommon_bound%var_ref)

        !! update the start solution
        ndim= (tcommon_bound%network%n_steps-1)*tcommon_bound%gap_saved_solutions+1
        do j=1,ndim-1
            tcommon_bound%start_solution(:,j)= tcommon_bound%start_solution(:,j+1)
        enddo 
        ! tcommon_bound%start_solution(1,ndim)= mass_flow
        ! tcommon_bound%start_solution(2:,ndim)= tcommon_bound%reduced_right_solution

        !! the input of the lstm unit is just the reduced order solution
        tcommon_bound%start_solution(:,ndim)= tcommon_bound%reduced_right_solution

        ! write(*,*) tcommon_bound%start_solution(:,ndim-1)
        ! write(*,*) tcommon_bound%start_solution(:,ndim)

    enddo



    ! call MPI_STOP


end subroutine

subroutine compute_common_boundary_value
    use mainvar
    use neuralnetworks
    implicit none
    type(common_bound),pointer:: tcommon_bound
    integer:: i,j

    !! evaluate the network

    do i=1,Ncommon_bound
        tcommon_bound => Common_bound_targets(i)

        !! get input 
        do j=1,tcommon_bound%network%n_steps
            tcommon_bound%network%input(:,j)= tcommon_bound%start_solution(:,(j-1)*tcommon_bound%gap_saved_solutions+1)
        enddo

        call forward_prop(tcommon_bound%network)

        ! write(*,*) tcommon_bound%network%output

        tcommon_bound%boundary_value= matmul(tcommon_bound%basis,tcommon_bound%network%output)

        tcommon_bound%boundary_value= tcommon_bound%boundary_value + tcommon_bound%var_ref

    enddo
    
    return
end subroutine
