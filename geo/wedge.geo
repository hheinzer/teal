n = 10;

Point(1) = {-2, -3, 0, 1.0};
Point(2) = {4, -3, 0, 1.0};
Point(3) = {4, 3, 0, 1.0};
Point(4) = {-2, 3, 0, 1.0};
Point(5) = {-1, 0, 0, 1.0};
Point(6) = {1, 0, 0, 1.0};
Point(7) = {0, 0.2, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 5};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7};
Plane Surface(1) = {1, 2};

Physical Curve("inflow", 8) = {1, 4};
Physical Curve("outflow", 9) = {2, 3};
Physical Curve("airfoil", 10) = {5, 6, 7};
Physical Surface("domain", 11) = {1};

Transfinite Curve {1, 2, 3, 4} = (n * 11 + 1) Using Progression 1;
Transfinite Curve {5} = (n * 12 + 1) Using Progression 1;
Transfinite Curve {6, 7} = (n * 7 + 1) Using Progression 1;
