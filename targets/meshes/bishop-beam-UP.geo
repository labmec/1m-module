SetFactory("OpenCASCADE");

el = 80;

L = 5;
a = 0.5;
b = 0.5;

Point(1) = {a, -b, 0, 1.0};
//+
Point(2) = {a, -b, L, 1.0};
//+
Point(3) = {a, b, L, 1.0};
//+
Point(4) = {a, b, 0, 1.0};
//+

Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+

//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Extrude {-2*a, 0, 0} {
  Surface{1}; 
}
//
Physical Volume("Domain", 1) = {1};
//+
Physical Surface("NormalDisplacement", 2) = {3};
//+
Physical Surface("NormalStress", 3) = {5,1,2,4,6};
//+
Physical Surface("TangentialDisplacement", 4) = {3};
//+
Physical Surface("TangentialStress", 5) = {5};
//+
Transfinite Curve {1,3,7,11} = el+1 Using Progression 1;
//+
Transfinite Curve {2,4,5,6,8,9,10,12} = Round(el/5)+1 Using Progression 1;
//+
Transfinite Surface {2} = {6, 4, 3, 5};
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Transfinite Surface {5} = {2, 8, 5, 3};
//+
Transfinite Surface {3} = {1, 4, 6, 7};
//+
Transfinite Surface {4} = {1, 7, 8, 2};
//+
Transfinite Surface {6} = {6, 5, 8, 7};
//+
Recombine Surface {2, 3, 6, 5, 1, 4};
//+
Transfinite Volume{1} = {1, 2, 3, 4, 7, 8, 5, 6};
//+
