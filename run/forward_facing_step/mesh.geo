SetFactory("OpenCASCADE");
nx = 600;
ny = 200;

Rectangle(1) = {0, 0, 0, 0.6, 0.2, 0};
Rectangle(2) = {0, 0.2, 0, 0.6, 0.8, 0};
Rectangle(3) = {0.6, 0.2, 0, 2.4, 0.8, 0};
Coherence;
Extrude {0, 0, 1} { Surface{1}; Surface{2}; Surface{3}; Layers {1}; Recombine; }

Transfinite Curve {1, 3, 6} = (nx / 3 * 0.6 + 1) Using Progression 1;
Transfinite Curve {8, 10} = (nx / 3 * 2.4 + 1) Using Progression 1;
Transfinite Curve {4, 2} = (ny / 1 * 0.2 + 1) Using Progression 1;
Transfinite Curve {7, 5, 9} = (ny / 1 * 0.8 + 1) Using Progression 1;
Transfinite Surface {1, 2, 3};
Recombine Surface {1, 2, 3};

Physical Surface("inlet", 29) = {7, 11};
Physical Surface("outlet", 30) = {14};
Physical Surface("symmetry", 31) = {15, 10};
Physical Surface("wall", 32) = {4, 5, 13};
Physical Volume("domain", 35) = {1, 2, 3};
