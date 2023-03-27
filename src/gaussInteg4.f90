!!======================================================
!!  Gauss integration functions
!!------------------------------------------------------
function GaussInteg_Tri( fun,xv1,xv2,xv3 )
    real*8 :: GaussInteg_Tri, xv1(2),xv2(2),xv3(2), &
          vx1(2),vx2(2),vx3(2),vx4(2),vx5(2),vx6(2),vx7(2),vx8(2),vx9(2),vx10(2),vx11(2),vx12(2),vx13(2)
    real*8,external :: fun
    
    x1=xv1(1); y1= xv1(2)
    x2=xv2(1); y2= xv2(2)
    x3=xv3(1); y3= xv3(2)
    area= 0.5*abs( (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1) )
    
    !! 7th order
    vx1= (xv1+xv2+xv3)/3.
    alpa= 0.479308067841920
    beta= 0.260345966079040
    vx2= alpa*xv1+ beta*(xv2+xv3)
    vx3= alpa*xv2+ beta*(xv1+xv3)
    vx4= alpa*xv3+ beta*(xv1+xv2)
    alpa= 0.869739794195568
    beta= 0.065130102902216
    vx5= alpa*xv1+ beta*(xv2+xv3)
    vx6= alpa*xv2+ beta*(xv1+xv3)
    vx7= alpa*xv3+ beta*(xv1+xv2)
    alpa= 0.048690315425316
    beta= 0.312865496004874
     vx8= alpa*xv1+ beta*xv2+(1-alpa-beta)*xv3
     vx9= alpa*xv1+ beta*xv3+(1-alpa-beta)*xv2
    vx10= alpa*xv2+ beta*xv1+(1-alpa-beta)*xv3
    vx11= alpa*xv2+ beta*xv3+(1-alpa-beta)*xv1
    vx12= alpa*xv3+ beta*xv1+(1-alpa-beta)*xv2
    vx13= alpa*xv3+ beta*xv2+(1-alpa-beta)*xv1
    
    GaussInteg_Tri= -0.149570044467682*fun(vx1)+ 0.175615257433208*(fun(vx2)+fun(vx3)+fun(vx4))+ &
                                    0.053347235608838*(fun(vx5)+ fun(vx6)+fun(vx7))+ &
                                    0.077113760890257*(fun(vx8)+fun(vx9)+fun(vx10)+fun(vx11)+fun(vx12)+fun(vx13))
    
    GaussInteg_Tri = GaussInteg_Tri* area
    ! write(*,*) GaussInteg_Tri, area, "Tri"
    return
end function

function GaussInteg_Quad(fun,xv1,xv2,xv3,xv4)
    real*8 :: GaussInteg_Quad,xv1(2),xv2(2),xv3(2),xv4(2)
    real*8,external :: fun,GaussInteg_Tri
    GaussInteg_Quad= GaussInteg_Tri(fun,xv1,xv2,xv3)+ &
                     GaussInteg_Tri(fun,xv1,xv3,xv4)
    return
end function
