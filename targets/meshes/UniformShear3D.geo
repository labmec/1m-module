// Gmsh project created on Fri Aug 25 11:33:46 2023
SetFactory("OpenCASCADE");

el = 2;
ndiv = el + 1;

L = 1.0;

Point(1) = {0.0, 0.0, 0.0, 1.0};
//+
Point(2) = {L, 0.0, 0.0, 1.0};
//+
Point(3) = {L, L, 0.0, 1.0};
//+
Point(4) = {0, L, 0.0, 1.0};
//+

Line(1) = {1, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {2, 1};
//+


//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, L} {
  Surface{1}; 
}

Physical Volume("Domain", 1) = {1};
//+
Physical Surface("ZeroNormalDisplacement", 2) = {3};
//+
Physical Surface("ZeroNormalStress", 3) = {1, 2, 4, 5, 6};
//+
Physical Surface("ZeroTangentialDisplacement", 4) = {3, 1, 6};
//+
Physical Surface("UnitTangentialStress", 5) = {5};
//+
Physical Surface("ZeroTangentialStress", 6) = {2, 4};

Transfinite Curve {4, 6, 7, 5, 3, 9, 12, 2, 1, 8, 11, 10} = ndiv Using Progression 1;
//+

//+
Transfinite Surface {2} = {2, 3, 5, 6};
//+
Transfinite Surface {5} = {4, 8, 5, 3};
//+
Transfinite Surface {4} = {4, 1, 7, 8};
//+
Transfinite Surface {1} = {1, 4, 3, 2};
//+
Transfinite Surface {6} = {8, 7, 6, 5};
//+
Transfinite Surface {3} = {1, 2, 6, 7};

Recombine Surface {1,2,3,4,5,6};
//+

Transfinite Volume {1};
