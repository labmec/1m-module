// Gmsh project created on Fri Aug 25 11:33:46 2023
SetFactory("OpenCASCADE");

elx = 1;
ely = 1;
elz = 1;
ndivx = elx + 1;
ndivy = ely + 1;
ndivz = elz + 1;

L = 1.0;

Point(1) = {0.0, 0.0, 0.0, 1.0};
//+
Point(2) = {L, 0.0, 0.0, 1.0};
//+
Point(3) = {L, L/2, 0.0, 1.0};
//+
Point(4) = {L, L, 0.0, 1.0};
//+
Point(5) = {0, L, 0.0, 1.0};
//+
Point(6) = {0, L/2, 0.0, 1.0};
//+

Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 6};
//+
Line(4) = {6, 1};
//+
Line(5) = {3, 4};
//+
Line(6) = {4, 5};
//+
Line(7) = {5, 6};
//+

Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {-3, 5, 6, 7};
//+
Plane Surface(2) = {2};
//+
Extrude {0, 0, -L} {
  Surface{1,2}; Layers{elz}; Recombine;
}
//+
Physical Volume("Domain1", 1) = {1};
//+
Physical Volume("Domain2", 2) = {2};
//+//+
Physical Surface("ZeroNormalDisp1", 3) = {1, 3, 4, 7};
//+
Physical Surface("ZeroNormalDisp2", 4) = {8, 11, 2};
//+
Physical Surface("UnitNormalDisp1", 5) = {5};
//+
Physical Surface("UnitNormalDisp2", 6) = {10};
//+
Physical Surface("ZeroNormalStress2", 7) = {9};
//+
Physical Surface("ZeroTangentialStress1", 8) = {1, 3, 7, 4, 5};
//+
Physical Surface("ZeroTangentialStress2", 9) = {9, 2, 11, 10, 8};
//+
Transfinite Curve {1, 12, 3, 15, 6, 19} = elx Using Progression 1;
//+
Transfinite Curve {7, 4, 10, 17, 5, 20, 14, 2} = ely Using Progression 1;
//+
Transfinite Curve {16, 18, 13, 8, 9, 11} = elz Using Progression 1;
//+
Transfinite Surface {3} = {6, 7, 8, 1};
//+
Transfinite Surface {8} = {6, 7, 11, 5}; 
//+
Transfinite Surface {5} = {2, 9, 10, 3};
//+
Transfinite Surface {10} = {3, 10, 12, 4};
//+
Transfinite Surface {9} = {4, 12, 11, 5};
//+
Transfinite Surface {4} = {1, 2, 9, 8};
//+
Transfinite Surface {1} = {1, 2, 3, 6};
//+
Transfinite Surface {2} = {6, 3, 4, 5};
//+
Transfinite Surface {11} = {7, 10, 12, 11};
//+
Transfinite Surface {7} = {8, 9, 10, 7};
//+
Transfinite Surface {6} = {6, 3, 10, 7};
//+
Recombine Surface {1,2,3,4,5,6,7,8,9,10,11};
//+
Transfinite Volume{2} = {6, 3, 4, 5, 11, 7, 10, 12};
//+
Transfinite Volume{1} = {1, 2, 3, 6, 7, 8, 9, 10};
