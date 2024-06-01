lc = 0.01;

Point(1) = {0, 0, 0, lc};
Point(2) = {0.5, 0, 0, lc};
Point(3) = {1, 0, 0, lc};
Point(4) = {1, 1, 0, lc};
Point(5) = {0.5, 1, 0, lc};
Point(6) = {0, 1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {2, 5};

Curve Loop(1) = {6, 1, 7, 5};
Plane Surface(1) = {1};
Curve Loop(2) = {7, -4, -3, -2};
Plane Surface(2) = {2};

Physical Curve("bottom", 8) = {1, 2};
Physical Curve("right", 9) = {3};
Physical Curve("top", 10) = {4, 5};
Physical Curve("left", 11) = {6};
Physical Surface("domain", 12) = {1, 2};

Recombine Surface {2};

