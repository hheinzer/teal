SetFactory("OpenCASCADE");
nx = 90;
ny = 30;
nz = 30;

Box(1) = {0, 0, 0, 9, 3, 3};

Periodic Surface {2} = {1} Translate {9, 0, 0};
Periodic Surface {6} = {5} Translate {0, 0, 3};

Transfinite Curve {11, 12, 9, 10} = (nx + 1) Using Progression 1;
Transfinite Curve {4, 2, 8, 6} = (ny + 1) Using Progression 1;
Transfinite Curve {1, 3, 7, 5} = (nz + 1) Using Progression 1;

Physical Volume("domain") = {1};
Physical Surface("left:right") = {1};
Physical Surface("right:left") = {2};
Physical Surface("bottom") = {3};
Physical Surface("top") = {4};
Physical Surface("back:front") = {5};
Physical Surface("front:back") = {6};
