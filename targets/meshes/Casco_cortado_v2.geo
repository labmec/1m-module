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
Mesh.MeshSizeMax = 15;
