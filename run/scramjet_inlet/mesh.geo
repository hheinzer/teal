SetFactory("OpenCASCADE");
lc = 0.02;

Point(1) = {0, 0, 0, lc};
Point(2) = {8, 0, 0, lc};
Point(3) = {8, 0.8, 0, lc};
Point(4) = {0, 1.7, 0, lc};
Point(5) = {2.1, 0.7, 0, lc};
Point(6) = {4.2, 0.22, 0, lc};
Point(7) = {6.7, 0.59, 0, lc};
Point(8) = {5.9, 0.7, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};
Extrude {0, 0, 1} { Surface{1}; Layers {1}; Recombine; }

Physical Surface("inlet", 25) = {5};
Physical Surface("outlet", 26) = {3};
Physical Surface("wall", 27) = {2, 4, 9, 6, 7, 8};
Physical Volume("domain", 28) = {1};
