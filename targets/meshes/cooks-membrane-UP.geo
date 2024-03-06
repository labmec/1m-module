SetFactory("OpenCASCADE");
a = 44;
b = 48;
c = 16;

scale = 0.25;

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
Physical Curve("ZeroNormalStress", 2) = {3, 1, 2};
//+
Physical Curve("ZeroNormalDisp", 3) = {4};
//+
Physical Curve("ZeroTangentialDisp", 4) = {4};
//+
Physical Curve("ShearStress", 5) = {2};
//+
//Transfinite Curve {1} = Round(scale*Hypot(a,b))+1 Using Progression 1;
//+
//Transfinite Curve {2} = Round(scale*c)+1 Using Progression 1;
//+
//Transfinite Curve {3} = Round(scale*Hypot(c,b))+1 Using Progression 1;
//+
//Transfinite Curve {4} = Round(scale*a)+1 Using Progression 1;

Transfinite Curve {1} = 1 Using Progression 1;
//+
Transfinite Curve {2} = 1 Using Progression 1;
//+
Transfinite Curve {3} = 1 Using Progression 1;
//+
Transfinite Curve {4} = 1 Using Progression 1;
//+
Transfinite Surface {1};
