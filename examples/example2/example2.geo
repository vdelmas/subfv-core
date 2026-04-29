//+
Point(1) = {-1, 0, 0, 1.0};
//+
Extrude {-2, 0, 0} {
  Point{1}; 
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Curve{1}; 
}
//+
Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/2} {
  Curve{2}; 
}
//+
Extrude {0, 0, 1} {
  Surface{5}; Surface{9}; Layers {1}; Recombine;
}
//+
Physical Surface("inflow_surf", 54) = {22, 44};
//+
Physical Surface("outflow_surf", 55) = {18, 48};
//+
Physical Volume("fluid", 56) = {1, 2};
//+
MeshSize {1, 2, 12, 11, 3, 21, 17, 4, 10, 28, 9, 32} = 0.02;
