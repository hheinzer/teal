SetFactory("OpenCASCADE");
Ds = 1.0;   // diameter of sphere
Rc = 5*Ds;  // radius of cylinder
Lu = 5*Ds;  // upstream extent of cylinder
Ld = 25*Ds; // downstream extent of cylinder

Sphere(1) = {0, 0, 0, Ds/2, -Pi/2, Pi/2, 2*Pi};
Cylinder(2) = {-Lu, 0, 0, Lu+Ld, 0, 0, Rc, 2*Pi};
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }

Point(1000) = {0, 0, 0};
Point(1001) = {Ld, 0, 0};
Line(1000) = {1000, 1001};

Field[1] = Distance;
Field[1].SurfacesList = {1};
Field[1].Sampling = 200;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].DistMin = 0.0;
Field[2].DistMax = Ds/3;
Field[2].Sigmoid = 1;
Field[2].SizeMin = 0.01;
Field[2].SizeMax = 0.2;

Field[3] = Distance;
Field[3].CurvesList = {1000};
Field[3].Sampling = 200;

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].DistMin = 0.0;
Field[4].DistMax = 7*Ds;
Field[4].Sigmoid = 1;
Field[4].SizeMin = 0.05;

Field[5] = Min;
Field[5].FieldsList = {2, 4};
Background Field = 5;

Physical Surface("wall", 7) = {1};
Physical Surface("farfield", 8) = {4, 2, 3};
Physical Volume("domain", 9) = {2};
