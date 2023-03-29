lc=0.0001;
L=1;
Point(1)={0,  0, 0,lc};
Point(2)={L,  0, 0,lc};
Point(3)={L,  lc,0,lc};
Point(4)={0  ,lc,0,lc};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Curve Loop(1)={1,2,3,4};
Plane Surface(1)={1};
Physical Curve(1)={1:4};
Physical Surface("My surface")={1:4};
//+
Transfinite Surface {1} Alternated;
//+
Transfinite Surface {1} Alternated;
