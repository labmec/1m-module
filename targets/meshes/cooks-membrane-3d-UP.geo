SetFactory("OpenCASCADE");
a = 44;
b = 48;
c = 16;

nx = 32;

Point(1) = {0.0, 0.0, 0.0, 1.0};
//+
Point(2) = {b, a, 0.0, 1.0};
//+
Point(3) = {b, a+c, 0.0, 1.0};
//+
Point(4) = {0.0, a, 0.0, 1.0};

Line(1) = {1,2};
//+
Line(2) = {2,3};
//+
Line(3) = {3,4};
//+
Line(4) = {4,1};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, -10} {
  Surface{1}; Layers{1}; Recombine;
}
//+

Physical Volume("Domain", 1) = {1};
//+
Physical Surface("ZeroNormalStress", 2) = {2,3,4};
//+
Physical Surface("ZeroNormalDisp", 3) = {1,5,6};
//+
Physical Surface("ZeroTangentialDisp", 4) = {5};
//+
Physical Surface("ShearStress", 5) = {3};
//+
//Physical Surface("ZeroShearStress", 6) = {1,2,4,6};
//+

Transfinite Curve {1} = nx+1 Using Progression 1;
//+
Transfinite Curve {2} = nx+1 Using Progression 1;
//+
Transfinite Curve {3} = nx+1 Using Progression 1;
//+
Transfinite Curve {4} = nx+1 Using Progression 1;
//+

Transfinite Surface {1} = {4, 1, 2, 3};
//+
Transfinite Surface {6} = {7, 6, 5, 8};
//+
Transfinite Surface {4} = {4, 3, 7, 8};
//+
Transfinite Surface {2} = {2, 1, 5, 6};
//+
Transfinite Surface {3} = {2, 6, 7, 3};
//+
Transfinite Surface {5} = {4, 8, 5, 1};
//+
Recombine Surface {2, 3, 6, 5, 1, 4};
//+
Transfinite Volume{1} = {2, 3, 4, 1, 5, 8, 7, 6};
//+
Physical Point("Sensor", 7) = {3};
//+