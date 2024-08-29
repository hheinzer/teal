lc = 0.02;

Point(1) = {0, 0, 0, lc};
Point(2) = {8, 0, 0, lc};
Point(3) = {8, 0.8, 0, lc};
Point(4) = {0, 2, 0, lc};
Point(5) = {2, 0.7, 0, lc};
Point(6) = {4, 0.2, 0, lc};
Point(7) = {7, 0.6, 0, lc};
Point(8) = {6, 0.7, 0, lc};

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
Extrude {0, 0, 1} {Surface{1}; Layers {1}; Recombine;}

Physical Surface("inlet", 51) = {33};
Physical Surface("outlet", 52) = {25};
Physical Surface("wall", 53) = {21, 29, 49, 45, 41, 37};
Physical Volume("domain", 54) = {1};
