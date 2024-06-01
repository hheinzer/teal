nx1 = 60;
nx2 = 240;
ny1 = 20;
ny2 = 80;

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {0.6, 0.0, 0.0, 1.0};
Point(3) = {0.6, 0.2, 0.0, 1.0};
Point(4) = {3.0, 0.2, 0.0, 1.0};
Point(5) = {3.0, 1.0, 0.0, 1.0};
Point(6) = {0.6, 1.0, 0.0, 1.0};
Point(7) = {0.0, 1.0, 0.0, 1.0};
Point(8) = {0.0, 0.2, 0.0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};
Line(9) = {3, 8};
Line(10) = {6, 3};

Curve Loop(1) = {1, 2, 9, 8};
Plane Surface(1) = {1};
Curve Loop(2) = {7, -9, -10, 6};
Plane Surface(2) = {2};
Curve Loop(3) = {5, 10, 3, 4};
Plane Surface(3) = {3};

Physical Curve("inflow", 11) = {8, 7};
Physical Curve("outflow", 12) = {4};
Physical Curve("symmetry", 13) = {1, 6, 5};
Physical Curve("wall", 14) = {2, 3};
Physical Surface("domain", 15) = {1, 2, 3};

Transfinite Curve {1, 9, 6} = nx1 + 1 Using Progression 1;
Transfinite Curve {3, 5} = nx2 + 1 Using Progression 1;
Transfinite Curve {8, 2} = ny1 + 1 Using Progression 1;
Transfinite Curve {7, 10, 4} = ny2 + 1 Using Progression 1;
Transfinite Surface {1, 2, 3};
Recombine Surface {1, 2, 3};
