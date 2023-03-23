!!======================================================
!!  Convseved variables to primitive variables
!!------------------------------------------------------
subroutine update_primitve_variables
    use mainvar
    implicit none
    integer:: i
    real*8 :: rho,u,v,p,t,Cp,gama,Rcpcv,Y(Nvar-Nflow+1)

    do i=1,NE

        rho =  WA(1,i)
        u   =  WA(2,i)/rho
        v   =  WA(3,i)/rho
        Y(1:Nvar-Nflow)=  WA(nflow+1:Nvar,i)/rho

        !! limiting the mass fraction
        ! Y(1:Nvar-Nflow)= max(0., Y(1:Nvar-Nflow))
        ! Y(1:Nvar-Nflow)= min(1., Y(1:Nvar-Nflow))

        !! throw the mass conservation error to the last species
        Y(Nvar-Nflow+1)= 1. - sum(Y(1:Nvar-Nflow))
        
        !! use Newton's method to compute temperature and pressure
        call Newton_temperature(WA(4,i)/rho-0.5*(u*u+v*v),Y,t,PA(4,i),Cp,gama,Rcpcv)
        p= rho*Rcpcv*t

        PA(1,i)= p
        PA(2,i)= u
        PA(3,i)= v
        PA(4,i)= t
        PA(nflow+1:Nvar,i)= Y(1:Nvar-Nflow)

        ! if(isnan(sum(PA(:,i)))) then
        !     write(*,*) 'error in update_primitve_variables'
        !     write(*,*) 'VARIABLES=',PA(:,i)
        !     write(*,*) 'WA=',WA(:,i)
        !     write(*,*) 'rho=',rho
        !     write(*,*) 'u=',u,'v=',v
        !     write(*,*) 'Y=',Y
        !     write(*,*) 'flux=',hWA(:,i)
        ! endif
    enddo

    call Connect_CellAverage

    return
end subroutine

subroutine primitve_to_conserved_variables
    use mainvar
    implicit none
    integer:: i
    real*8 :: pp,uu,vv,tt,rr,YY(Nvar-Nflow+1),Cp,gama,Rcpcv,entha(Nvar-Nflow+1),hh

    do i=1,NE
      pp=PA(1,i)
      uu=PA(2,i)
      vv=PA(3,i)
      tt=PA(4,i)
      YY(1:Nvar-Nflow)=PA(nflow+1:NVar,i); YY(Nvar-Nflow+1)= 1. - sum(YY(1:Nvar-Nflow))

      ! if(tt<200 .or. tt>5000) then
      !   write(*,*) 'error in temperature, tt=', tt, 'of cell i=',i
      !   write(*,*) 'p,u,v,t,YY=',pp,uu,vv,tt,YY(:)
      !   stop
      ! endif
      
      call ComputeGasParameter(YY, tt, Cp, gama, Rcpcv)
      rr=pp/(Rcpcv*tt)

      call compute_enthalpy(tt,entha)
      hh= sum(YY(:)*entha(:))

      WA(1,i)= rr
      WA(2,i)= rr*uu
      WA(3,i)= rr*vv
      WA(4,i)= rr*(hh + 0.5*(uu**2+vv**2) ) - pp
      WA(nflow+1:NVar,i)= rr*YY(1:Nvar-Nflow)

    enddo

    return
end subroutine


!!================================================================
!!  Update ghost cell variables based on primitive variables
!!----------------------------------------------------------------
subroutine update_ghostcell_variables
    use mainvar
    use neuralnetworks
    implicit none
    real*8 :: rl,ul,vl,tl,pl,rr,ur,vr,tr,pr,&
              sav1,sav2,sav,sav1n,sav2n,vn,&
              x0,flag,xmid,ymid,ac,Cp,gama,Rcpcv,Yl(Nvar-Nflow+1),Yr(Nvar-Nflow+1)
    real*8 :: xv1(2),xv2(2),xv3(2),xv4(2),xy(2)
    integer:: i,j,iface,id_face,id_face_number
    type(common_bound),pointer:: tcommon_bound

    do j=1,Nghostcell

        i= Hostcell(j)

        pl= PA(1,i)
        ul= PA(2,i)
        vl= PA(3,i)
        tl= PA(4,i)
        Yl(1:Nvar-Nflow)=PA(Nflow+1:Nvar,i); Yl(Nvar-Nflow+1)= 1. - sum(Yl(1:Nvar-Nflow))
        call ComputeGasParameter(Yl, tl, Cp, gama, Rcpcv)
        rl= pl/(Rcpcv*tl)


        xv1=Coor_Ghostnode(:,1,j)
        xv2=Coor_Ghostnode(:,2,j)
        xv3=Coor_Ghostnode(:,3,j)
        xv4=Coor_Ghostnode(:,4,j)
        
        xmid=0.5*(xv1(1)+xv2(1))
        ymid=0.5*(xv1(2)+xv2(2))
        xy=(/xmid,ymid/)

        sav1=   xv2(2)-xv1(2)  
        sav2= -(xv2(1)-xv1(1))

        sav = max(sqrt(sav1*sav1+sav2*sav2),small)
        sav1n=sav1/sav
        sav2n=sav2/sav

        ac= sqrt(gama*pl/rl)
        
        select case(BCtype_Ghostcell(j))
 
        case(1) !! viscous solid wall
            
            pr=pl; tr= tl; Yr= Yl
            ur=-ul;vr=-vl

        case(2) !! symmetry

            pr= pl; tr=tl; Yr= Yl
            vn= ul*sav1n + vl*sav2n
            ur= ul-2*vn*sav1n
            vr= vl-2*vn*sav2n

        case(3) !! oxidizer inlet

            call oxidizer_inlet_bc(pl,pr,ur,vr,tr,rr,Yr,Cp,gama,Rcpcv)

            ! pr= 2*pr - pl
            ! ur= 2*ur - ul
            ! vr= 2*vr - vl
            ! tr= 2*tr - tl
            ! Yr= 2*Yr - Yr

        case(4) !! fuel inlet
        
            call fuel_inlet_bc(pl,pr,ur,vr,tr,rr,Yr,Cp,gama,Rcpcv)

            ! pr= 2*pr - pl
            ! ur= 2*ur - ul
            ! vr= 2*vr - vl
            ! tr= 2*tr - tl
            ! Yr= 2*Yr - Yl

        case(5) !! nozzle outlet

            vn= ul*sav1n + vl*sav2n

            if(abs(vn)>ac) then  !! supersonic
                pr= pl
                ur= ul
                vr= vl
                tr= tl
                Yr= Yl
            else                  !! subsonic
                ! pr= p_out
                ! ur= ul + sav1n*(pl-pr)/(rl*ac)
                ! vr= vl + sav2n*(pl-pr)/(rl*ac)
                ! rr= rl + (pr-pl)/(ac**2)
                ! tr= pr/(Rcpcv*rr)

                ! pr= 2*pr - pl
                ! ur= 2*ur - ul
                ! vr= 2*vr - vl
                ! tr= 2*tr - tl
                ! Yr= Yl

                pr= p_out
                ur= ul
                vr= vl
                tr= tl
                Yr= Yl
            endif


            !! supersonic outlet
            ! pr= pl
            ! ur= ul
            ! vr= vl
            ! tr= tl
            ! Yr= Yl

            !! static pressure 2
            ! pr= p_out
            ! ur= ul + sav1n*(pl-pr)/(rl*ac)
            ! vr= vl + sav2n*(pl-pr)/(rl*ac)
            ! rr= rl + (pr-pl)/(ac**2)
            ! tr= pr/(Rcpcv*rr)

            ! pr= 2*pr - pl
            ! ur= 2*ur - ul
            ! vr= 2*vr - vl
            ! tr= 2*tr - tl

            ! Yr= Yl

        case(6) !! inviscid solid wall

            pr= pl; tr=tl; Yr= Yl
            vn= ul*sav1n + vl*sav2n
            ur= ul-2*vn*sav1n
            vr= vl-2*vn*sav2n

        case(1000:) !!  common boundary
            
            iface= Hostface(j)

            id_face= index_common_boundary(1,iface)
            id_face_number= index_common_boundary(2,iface)

            if(id_face==0 .or. id_face_number==0) then
              write(*,*) 'error in finding the common boundary face'
              call MPI_STOP
            endif

            tcommon_bound => Common_bound_targets(id_face) 

            pr= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+1)
            ur= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+2)
            vr= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+3)
            tr= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+4)

            Yr(1:Nvar-Nflow)= tcommon_bound%boundary_value(Nvar*(id_face_number-1)+Nflow+1:Nvar*id_face_number)
            Yr(Nvar-Nflow+1)= 1. - sum(Yr(1:Nvar-Nflow))

            ! call ComputeGasParameter(Yr,tr,Cp,gama,Rcpcv)
            ! rr= pr/(Rcpcv*tr)

            pr= 2*pr - pl
            ur= 2*ur - ul
            vr= 2*vr - vl
            tr= 2*tr - tl
            Yr= 2*Yr - Yl
            
        case default
            write(*,*) 'error in BCtype of ghost cess'
            write(*,*) 'bc type is:',BCtype_Ghostcell(j)
            call MPI_STOP
        end select

        GPA(1,j)= pr
        GPA(2,j)= ur
        GPA(3,j)= vr
        GPA(4,j)= tr
        GPA(nflow+1:Nvar,j)= Yr(1:Nvar-Nflow)

    enddo
end subroutine