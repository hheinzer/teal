SetFactory("OpenCASCADE");
n = 50;

Box(1) = {-1, -1, -1, 2, 2, 2};

Periodic Surface {2} = {1} Translate {2, 0, 0};
Periodic Surface {4} = {3} Translate {0, 2, 0};
Periodic Surface {6} = {5} Translate {0, 0, 2};

Physical Surface("left", 13) = {1};
Physical Surface("right", 14) = {2};
Physical Surface("bottom", 15) = {3};
Physical Surface("top", 16) = {4};
Physical Surface("back", 17) = {5};
Physical Surface("front", 18) = {6};
Physical Volume("domain", 19) = {1};

Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = (n + 1) Using Progression 1;

Mesh.Algorithm3D = 10; // HXT
Mesh.OptimizeNetgen = 1; // improve quality of tetrahedral elements
General.NumThreads = 0; // use system default, i.e. OMP_NUM_THREADS
