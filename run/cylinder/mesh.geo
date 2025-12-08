SetFactory("OpenCASCADE");
inner = 1;
outer = 10;
length = 2;
n = 100;

Rectangle(1) = {0, -outer / 2, 0, length * outer, outer, 0};
Disk(2) = {0, 0, 0, inner / 2, inner / 2};
Disk(3) = {0, 0, 0, outer / 2, outer / 2};
BooleanUnion{ Surface{1}; Delete; }{ Surface{3}; Delete; }
BooleanDifference{ Surface{1}; Delete; }{ Surface{2}; Delete; }
Extrude {0, 0, 1} { Surface{1}; Layers {1}; Recombine; }

Physical Surface("wall", 20) = {6};
Physical Surface("farfield", 21) = {4, 2, 3, 5};
Physical Volume("domain", 22) = {1};

Transfinite Curve {5} = 2 * (n + 1) Using Progression 1;
Transfinite Curve {8, 7} = (n + 1) / 2 Using Progression 1;
Transfinite Curve {6, 9} = length * (n + 1) / 2 Using Progression 1;
