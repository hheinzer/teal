SetFactory("OpenCASCADE");
nx = 180;
ny = 60;
nz = 20;

Box(1) = {0, 0, 0, 9, 3, 1};

Periodic Surface {2} = {1} Translate {9, 0, 0};

Physical Volume("domain", 13) = {1};
Physical Surface("left:right", 14) = {1};
Physical Surface("right:left", 15) = {2};
Physical Surface("bottom", 16) = {3};
Physical Surface("top", 17) = {4};

Transfinite Curve {11, 12, 9, 10} = (nx + 1) Using Progression 1;
Transfinite Curve {4, 2, 8, 6} = (ny + 1) Using Progression 1;
Transfinite Curve {1, 3, 7, 5} = (nz + 1) Using Progression 1;
