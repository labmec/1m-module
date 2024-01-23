SetFactory("OpenCASCADE");

el = 10;

L = 5;
a = 0.5;
b = 0.5;

Point(1) = {-a, -b, 0, 1.0};
//+
Point(2) = {-a, -b, L, 1.0};
//+
Point(3) = {-a, b, L, 1.0};
//+
Point(4) = {-a, b, 0, 1.0};
//+

Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+

Curve Loop(1) = {4,1,2,3};
//+
Plane Surface(1) = {1};
//+
Extrude {2*a, 0, 0} {
  Surface{1}; Layers{el/5}; Recombine;
}
//+

Physical Volume("Domain", 1) = {1};
//+
Physical Surface("ZeroNormalDisplacement", 2) = {2};
//+
Physical Surface("EndForce", 3) = {4};
//+
Physical Surface("ZeroNormalStress", 4) = {1,3,5,6};
//+
Physical Surface("ZeroTangentialDisplacement", 5) = {2};
//+
Physical Surface("EndShear", 6) = {4};
//+
Physical Surface("ZeroTangentialStress", 7) = {1,3,5,6};
//+
Transfinite Curve {1,3,9,12} = el+1 Using Progression 1;
//+
Transfinite Curve {2,4,5,6,7,8,10,11} = Round(el/5)+1 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 3, 4} Left;
//+
Transfinite Surface {2} = {4, 5, 6, 1} Left;
//+
Transfinite Surface {3} = {1, 6, 7, 2} Left;
//+
Transfinite Surface {4} = {3, 2, 7, 8} Left;
//+
Transfinite Surface {5} = {5, 4, 3, 8} Left;
//+
Transfinite Surface {6} = {5, 8, 7, 6} Right;
//+
Recombine Surface {2, 3, 6, 5, 1, 4};
//+
Transfinite Volume{1} = {1, 2, 3, 4, 6, 7, 8, 5};
//+

