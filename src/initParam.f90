
subroutine InitParam
    use mainvar
    implicit none
    integer :: i,nd,ig1,ig2
    real*8  :: psi,alpha

!!-----------------------------
!! face integration
!!-----------------------------    
!! 2nd
    wgg = (/ 1.  /)                             !! function value weighted
    wxx = (/ 0.5 /)  

!!-------------------------------
!! volume integration
!!-------------------------------
    !!NG=1, triangle, 2nd order
    wxxV(:,1)= (/ 1./3, 1./3, 1./3 /)
    wggV(:)  = (/ 1. /)
    
    !! NG_V4=1, quadrilateral, 2nd-order
    wxxV4(:,1)= (/ 0.5, 0.5 /)
    wggV4(:)  = (/ 1. /)


!!-------------------------------
!! time integration
!!-------------------------------

    if(IfSteady) then
        write(*,*) 'We are solving a steady-state problem'
    else
        write(*,*) 'We are solving a unsteady problem'
    endif

    !!--------------------------------
    !! Explicit Methods
    !!--------------------------------
    if( IfExplicit )then
        select case(InLoop)
        case(1)     !! Forward Euler
            cs(1)=0.
        case(3)     !! explicit Runge-Kutta 3
            cs(1)=0.
            cs(2)=1.
            cs(3)=0.5
        case default
            write(*,*) 'Error: No such explicit Rnge-Kutta time-marching setting.'
            write(*,*) 'InLoop=',InLoop
            call MPI_STOP
        end select

        InnerLoop=1

        write(*,*) 'using an explicit Runge-Kutta method, with ',InLoop, ' stages', InnerLoop,' inner iterations'

    !!--------------------------------
    !! Implicit Methods
    !!--------------------------------
    else

        !! steady problems
        if( IfSteady )then
            InLoop= 1
            InnerLoop=1

            write(*,*) 'using a LU-SGS method, with ',InLoop, ' stages', InnerLoop,' inner iterations'

        !! unsteady problems
        else              
            dt= totaltime/nstep_max 

            select case(implicit_time_method)

            case(1)   !! implicit Runge-Kutta p, s stages

                ! p=2, s=InLoop=1  !p=4, s=InLoop=3 

                !! implicit RK coefficients
                asj= 0.
                bs = 0.
                cs = 0.

                select case(InLoop)
                case(1)
                    asj(1,1)= 0.5
                    bs(1)= 1.0

                case(3)
                    psi=0.128886400515720
                    alpha=2.*psi-1
                    asj(1,1)= psi 
                    asj(2,1:2)= (/ 0.5-psi,psi /)
                    asj(3,1:3)= (/ 2.*psi,1.-4.*psi,psi /)

                    bs(1)=1./6./alpha**2
                    bs(2)=1.-2.*bs(1)
                    bs(3)=bs(1)
                case default
                    write(*,*) 'Error: No such implicit Runge-Kutta time-marching setting.'
                    write(*,*) 'InLoop=',InLoop
                    call MPI_STOP
                end select


                do ig1=1,InLoop
                    cs(ig1)= sum(asj(ig1,1:InLoop))
                enddo

                write(*,*) 'using an implicit Runge-Kutta method, with ',InLoop, ' stages', InnerLoop,' inner iterations'
            case(2)     !! Trapezoidal Rule
                InLoop= 1
                cs(1)= 1.

                write(*,*) 'using an implicit Trapezoidal Rule method, with ',InLoop, ' stages', InnerLoop,' inner iterations'
            case(3)     !! Implicit Euler
                InLoop= 1
                cs(1)= 1.

                write(*,*) 'using an implicit Euler method, with ',InLoop, ' stages', InnerLoop,' inner iterations'
            case(4)     !! Implicit Midpoint
                InLoop= 1
                cs(1)= 1.

                write(*,*) 'using an implicit midpoint method, with ',InLoop, ' stages', InnerLoop,' inner iterations'
            case(5)     !! BDF2
                InLoop= 1
                cs(1)= 1.

                write(*,*) 'using an BDF2 method, with ',InLoop, ' stages', InnerLoop,' inner iterations'

            case default
                write(*,*) 'Error: No such implicit time stepping method.'
                write(*,*) 'implicit_time_method=',implicit_time_method
                stop
            end select

            write(*,*) 'with a time step size of dt=', dt
        endif 
        
    endif


    
    return
end subroutine

