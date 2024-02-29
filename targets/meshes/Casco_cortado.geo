SetFactory("OpenCASCADE");
//+
Merge "Casco_cortado.igs";
//+
Physical Volume("dom", 1) = {1};
//+
Physical Surface("pressure", 2) = {15,150};
//+
Physical Surface("sim_x", 3) = {45};
//+
Physical Surface("sim_y", 4) = {208};
//+
Physical Surface("top", 5) = {171};
//+
Physical Surface("bottom", 6) = {152};
//+
Mesh.MeshSizeMax = 15;
