SetFactory("OpenCASCADE");
nx = 300;
ny = 100;

Rectangle(1) = {0, 0, 0, 3, 1, 0};
Rectangle(2) = {0.6, 0, 0, 2.4, 0.2, 0};
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
Extrude {0, 0, 1} { Surface{1}; Layers {1}; Recombine; }

Physical Surface("inlet", 25) = {2};
Physical Surface("outlet", 26) = {4};
Physical Surface("symmetry", 27) = {3};
Physical Surface("wall", 28) = {7, 6, 5};
Physical Volume("domain", 31) = {1};

Transfinite Curve {10} = (nx + 1) Using Progression 1;
Transfinite Curve {9} = (ny + 1) Using Progression 1;
Transfinite Curve {12} = (nx / 3 * 0.6 + 1) Using Progression 1;
Transfinite Curve {7} = (nx / 3 * 2.4 + 1) Using Progression 1;
Transfinite Curve {8} = (ny / 1 * 0.2 + 1) Using Progression 1;
Transfinite Curve {11} = (ny / 1 * 0.8 + 1) Using Progression 1;
