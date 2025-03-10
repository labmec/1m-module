
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
Physical Surface("Domain", 1) = {1};
//+
Physical Curve("FixedEnd", 2) = {4};
//+
Physical Curve("ShearStress", 3) = {2};
//+
Physical Point("Sensor", 4) = {3};
//+

Transfinite Curve {1} = nx+1 Using Progression 1;
//+
Transfinite Curve {2} = nx+1 Using Progression 1;
//+
Transfinite Curve {3} = nx+1 Using Progression 1;
//+
Transfinite Curve {4} = nx+1 Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
