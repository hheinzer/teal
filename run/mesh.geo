SetFactory("OpenCASCADE");

Box(1) = {0, 0, -0.5, 9, 3, 1};

Physical Volume("domain") = {1};
Physical Surface("bottom") = {3};
Physical Surface("top") = {4};
Periodic Surface {2} = {1} Translate {9, 0, 0};
Periodic Surface {6} = {5} Translate {0, 0, 1};
