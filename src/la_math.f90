subroutine pseudo_inverse_lapack(a,m,n,b)
    implicit none
    external:: dgesvd
    integer, intent(in):: m,n
    !    .. Parameters ..
    Integer, Parameter :: nb = 64
    !     .. Local Scalars ..
    Integer :: i, info, lda, ldu, ldvt, lwork, mndim
    !     .. Local Arrays ..
    double precision:: a(m,n),b(n,m)
    ! double precision:: acopy(m,n)
    double precision, Allocatable :: s(:), u(:, :), vt(:, :), work(:)
    double precision :: dummy(1)

    ! acopy= a

    mndim= min(m,n)

    lda = m
    ldu = m
    ldvt = n
    Allocate (s(min(m,n)), vt(mndim,n), u(m,mndim));s=0.;vt=0.;u=0.

    !     Use routine workspace query to get optimal workspace.
    lwork = -1
    Call dgesvd('S', 'S', m, n, a, lda, s, u, ldu, vt, ldvt, dummy, lwork, &
    info)

    !     Make sure that there is enough workspace for block size nb.
    ! lwork = max(m+4*n+nb*(m+n), nint(dummy(1)))
    lwork = nint(dummy(1))
    Allocate (work(lwork));work=0.

    !     Compute the singular values and left and right singular vectors
    !     of A.

    Call dgesvd('S', 'S', m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, &
    info)

    if(info/=0) then
        write(*,*) 'error in routine dgesvd'
        write(*,*) 'error type, info=',info
        write(*,*) 'm=',m,'n=',n
        write(*,*) 's=',s
        write(*,*) 'lwork=',lwork

        ! write(*,*) 'a1=',acopy(1,:)
        ! write(*,*) 'a2=',acopy(2,:)
        ! write(*,*) 'a3=',acopy(3,:)
        stop
    endif

    do i=1,mndim
        vt(i,:)= vt(i,:)/s(i)
    enddo

    b= matmul(transpose(vt),transpose(u))


    deallocate(work,s,u,vt)

end subroutine

subroutine matrix_inverse_lapack(a,n)
    implicit none
    external:: dgetrf,dgesvd
    !     .. Local Scalars ..
    integer, intent(in):: n
    integer :: info, lda, lwork, nb
    !     .. Local Arrays ..
    double precision:: a(n, n)
    double precision,allocatable:: work(:)
    double precision:: dummy(1)
    integer:: ipiv(n)

    lda = n

    !     Factorize A

    Call dgetrf(n, n, a, lda, ipiv, info)
    if(info/=0) then
        write(*,*) 'error in routine dgetrf'; stop
    endif

    !     Compute inverse of A
    lwork=-1
    Call dgetri(n, a, lda, ipiv, dummy, lwork, info)

    lwork= nint(dummy(1))
    Allocate (work(lwork));work=0.

    Call dgetri(n, a, lda, ipiv, work, lwork, info)
    if(info/=0) then
        write(*,*) 'error in routine dgetri'; stop
    endif

    deallocate(work)

end subroutine

subroutine linear_solve_lapack(a,n,b)
    implicit none
    external:: dgesv
    !     .. Local Scalars ..
    integer, intent(in):: n
    integer :: info
    !     .. Local Arrays ..
    double precision:: a(n, n),b(n),x(n)
    integer:: ipiv(n)

    !     Solve the equations Ax = b for x

    Call dgesv(n, 1, a, n, ipiv, b, n, info)
    if(info/=0) then
        write(*,*) 'error in routine dgesv'; stop
    endif

end subroutine

subroutine advanced_linear_solve_lapack(a,n,b,x)
    implicit none
    external:: dgesvx
    !     .. Local Scalars ..
    integer, intent(in):: n
    integer :: info
    !     .. Local Arrays ..
    double precision:: a(n,n),a_copy(n,n), af(n,n), b(n), berr, c(n), &
        ferr, r(n), work(4*n), x(n),rcond
    integer:: ipiv(n),iwork(n)
    Character (1) :: equed

    a_copy= a

!     Solve the equations AX = B for X

    Call dgesvx('E', 'N', n, 1, a, n, af, n, &
    ipiv, equed, r, c, b, n, x, n, rcond, ferr, berr, work, iwork, &
    info)

    ! If ((info==0) .Or. (info==n+1)) Then

    !     Write (*, *)
    !     Write (*, *) 'Backward errors (machine-dependent)'
    !     Write (*, *) berr
    !     Write (*, *)
    !     Write (*, *) 'Estimated forward error bounds (machine-dependent)'
    !     Write (*, *) ferr
    !     Write (*, *)
    !     If (equed=='N') Then
    !       Write (*, *) 'A has not been equilibrated'
    !     Else If (equed=='R') Then
    !       Write (*, *) 'A has been row scaled as diag(R)*A'
    !     Else If (equed=='C') Then
    !       Write (*, *) 'A has been column scaled as A*diag(C)'
    !     Else If (equed=='B') Then
    !       Write (*, *) &
    !         'A has been row and column scaled as diag(R)*A*diag(C)'
    !     End If
    !     Write (*, *)
    !     Write (*, *) &
    !       'Reciprocal condition number estimate of scaled matrix'
    !     Write (*, *) rcond
    !     Write (*, *)
    !     Write (*, *) 'Estimate of reciprocal pivot growth factor'
    !     Write (*, *) work(1)

    !     If (info==n+1) Then
    !       Write (*, *)
    !       Write (*, *) 'The matrix A is singular to working precision'
    !     End If
    ! Else
    !     Write (*, *) 'The (', info, ',', info, ')', &
    !   ' element of the factor U is zero'
    ! End If

    if(info/=0) then
        write(*,*) 'error in routine dgesvx'

        Write (*, *)
        Write (*, *) 'Backward errors (machine-dependent)'
        Write (*, *) berr
        Write (*, *)
        Write (*, *) 'Estimated forward error bounds (machine-dependent)'
        Write (*, *) ferr
        Write (*, *)
        If (equed=='N') Then
          Write (*, *) 'A has not been equilibrated'
        Else If (equed=='R') Then
          Write (*, *) 'A has been row scaled as diag(R)*A'
        Else If (equed=='C') Then
          Write (*, *) 'A has been column scaled as A*diag(C)'
        Else If (equed=='B') Then
          Write (*, *) &
            'A has been row and column scaled as diag(R)*A*diag(C)'
        End If
        Write (*, *)
        Write (*, *) &
          'Reciprocal condition number estimate of scaled matrix'
        Write (*, *) rcond
        Write (*, *)
        Write (*, *) 'Estimate of reciprocal pivot growth factor'
        Write (*, *) work(1)

        If (info==n+1) Then
          Write (*, *)
          Write (*, *) 'The matrix A is singular to working precision'
        End If

        Write(*,*) 'The matrix is', a_copy
        stop
    endif

end subroutine