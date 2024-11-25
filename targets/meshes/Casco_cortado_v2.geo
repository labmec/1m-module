SetFactory("OpenCASCADE");

Merge "Casco_cortado_v2.igs";
//+
Physical Volume("dom", 1) = {1};
//+
Physical Surface("pressure", 2) = {68,23};
//+
Physical Surface("sim_x", 3) = {56};
//+
Physical Surface("sim_y", 4) = {54};
//+
Physical Surface("top", 5) = {67};
//+
Physical Surface("bottom", 6) = {53};
//+
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
//+
Field[1] = Threshold;
//+
Delete Field [1];
//+
Field[1] = Box;
//+
Field[1].VIn = 1000;
//+
Field[1].VOut = 20;
//+
Field[1].XMax = 1100;
//+
Field[1].XMin = -1100;
//+
Field[1].YMax = 1100;
//+
Field[1].YMin = -1100;
//+
Field[1].ZMax = 1100;
//+
Field[1].ZMin = -1100;
//+
Background Field = 1;
//+
Physical Curve("stiffner", 7) = {118, 138, 134, 130, 126};
