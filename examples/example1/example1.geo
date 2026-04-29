//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Physical Surface("left_surf", 13) = {1};
//+
Physical Surface("right_surf", 14) = {2};
//+
Physical Volume("fluid", 15) = {1};
