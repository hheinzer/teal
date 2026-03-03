SetFactory("OpenCASCADE");
n = 6;

Box(1) = {0, 0, 0, 9, 3, 3};

Physical Surface("left") = {1};
Physical Surface("right") = {2};
Physical Surface("bottom") = {3};
Physical Surface("top") = {4};
Physical Surface("back") = {5};
Physical Surface("front") = {6};
Physical Volume("domain") = {1};

Periodic Surface {2} = {1} Translate {9, 0, 0};
Periodic Surface {6} = {5} Translate {0, 0, 3};

Transfinite Curve {9, 10, 12, 11} = ((3 * n) + 1) Using Progression 1;
Transfinite Curve {4, 2, 8, 6} = (n + 1) Using Progression 1;
Transfinite Curve {7, 5, 3, 1} = (n + 1) Using Progression 1;
