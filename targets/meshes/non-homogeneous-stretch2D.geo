// Gmsh project created on Fri Aug 25 11:33:46 2023
SetFactory("OpenCASCADE");

elx = 1;
ely = 1;
ndivx = elx + 1;
ndivy = ely + 1;

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
Physical Surface("Domain1", 1) = {1};
//+

Curve Loop(2) = {-3, 5, 6, 7};
//+
Plane Surface(2) = {2};
//+
Physical Surface("Domain2", 2) = {2};

Physical Curve("ZeroNormalDisp1", 3) = {1,4};
//+
Physical Curve("ZeroNormalDisp2", 4) = {7};
//+
Physical Curve("UnitNormalDisp1", 5) = {2};
//+
Physical Curve("UnitNormalDisp2", 6) = {5};
//+
Physical Curve("ZeroNormalStress2", 7) = {6};
//+
Physical Curve("ZeroTangentialStress1", 8) = {1,2,4};
//+
Physical Curve("ZeroTangentialStress2", 9) = {5,6,7};
//+

Transfinite Surface {1} = {1, 2, 3, 6};
//+
Transfinite Surface {2} = {6, 3, 4, 5};
//+
Transfinite Curve {1} = ndivx Using Progression 1;
//+
Transfinite Curve {2} = ndivy Using Progression 1;
//+
Transfinite Curve {3} = ndivx Using Progression 1;
//+
Transfinite Curve {4} = ndivy Using Progression 1;
//+
Transfinite Curve {5} = ndivy Using Progression 1;
//+
Transfinite Curve {6} = ndivx Using Progression 1;
//+
Transfinite Curve {7} = ndivy Using Progression 1;
//+

Recombine Surface {1,2};
