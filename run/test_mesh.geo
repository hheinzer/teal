SetFactory("OpenCASCADE");
nx = 3;
ny = 3;

Rectangle(1) = {0, 0, 0, 9, 3, 0};
Extrude {0, 0, 1} {
  Surface{1}; Layers {1}; Recombine;
}

Physical Volume("domain", 13) = {1};
Physical Surface("left", 14) = {5};
Physical Surface("right", 15) = {3};
Physical Surface("bottom", 16) = {2};
Physical Surface("top", 17) = {4};

Transfinite Curve {1, 3} = (nx + 1) Using Progression 1;
Transfinite Curve {2, 4} = (ny + 1) Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};
