SetFactory("OpenCASCADE");
a = 30; // sweep angle
b = 1.47601797621980; // semispan
c = 0.562; // tip chord
r = 10; // sphere radius

Point(newp) = { 0.0000000, 0.0000000, 0};
Point(newp) = { 0.0000165, 0.0006914, 0};
Point(newp) = { 0.0000696, 0.0014416, 0};
Point(newp) = { 0.0001675, 0.0022554, 0};
Point(newp) = { 0.0003232, 0.0031382, 0};
Point(newp) = { 0.0005508, 0.0040959, 0};
Point(newp) = { 0.0008657, 0.0051343, 0};
Point(newp) = { 0.0012868, 0.0062598, 0};
Point(newp) = { 0.0018364, 0.0074784, 0};
Point(newp) = { 0.0025441, 0.0087958, 0};
Point(newp) = { 0.0034428, 0.0102163, 0};
Point(newp) = { 0.0045704, 0.0117419, 0};
Point(newp) = { 0.0059751, 0.0133708, 0};
Point(newp) = { 0.0077112, 0.0150951, 0};
Point(newp) = { 0.0098413, 0.0168984, 0};
Point(newp) = { 0.0124479, 0.0187537, 0};
Point(newp) = { 0.0156171, 0.0206220, 0};
Point(newp) = { 0.0194609, 0.0224545, 0};
Point(newp) = { 0.0241067, 0.0242004, 0};
Point(newp) = { 0.0297008, 0.0258245, 0};
Point(newp) = { 0.0364261, 0.0273317, 0};
Point(newp) = { 0.0444852, 0.0287912, 0};
Point(newp) = { 0.0541248, 0.0303278, 0};
Point(newp) = { 0.0656303, 0.0320138, 0};
Point(newp) = { 0.0793366, 0.0338372, 0};
Point(newp) = { 0.0956354, 0.0357742, 0};
Point(newp) = { 0.1149796, 0.0377923, 0};
Point(newp) = { 0.1378963, 0.0398522, 0};
Point(newp) = { 0.1649976, 0.0419089, 0};
Point(newp) = { 0.1919327, 0.0436214, 0};
Point(newp) = { 0.2187096, 0.0450507, 0};
Point(newp) = { 0.2453310, 0.0462358, 0};
Point(newp) = { 0.2717978, 0.0471987, 0};
Point(newp) = { 0.2981113, 0.0479494, 0};
Point(newp) = { 0.3242726, 0.0484902, 0};
Point(newp) = { 0.3502830, 0.0488183, 0};
Point(newp) = { 0.3761446, 0.0489296, 0};
Point(newp) = { 0.4018567, 0.0488202, 0};
Point(newp) = { 0.4274223, 0.0484833, 0};
Point(newp) = { 0.4528441, 0.0479351, 0};
Point(newp) = { 0.4781197, 0.0471661, 0};
Point(newp) = { 0.5032514, 0.0461903, 0};
Point(newp) = { 0.5282426, 0.0450209, 0};
Point(newp) = { 0.5530937, 0.0436741, 0};
Point(newp) = { 0.5778043, 0.0421684, 0};
Point(newp) = { 0.6023757, 0.0405241, 0};
Point(newp) = { 0.6268104, 0.0387613, 0};
Point(newp) = { 0.6511093, 0.0368990, 0};
Point(newp) = { 0.6752726, 0.0349542, 0};
Point(newp) = { 0.6993027, 0.0329402, 0};
Point(newp) = { 0.7231995, 0.0308662, 0};
Point(newp) = { 0.7469658, 0.0287365, 0};
Point(newp) = { 0.7705998, 0.0265505, 0};
Point(newp) = { 0.7941055, 0.0243027, 0};
Point(newp) = { 0.8174828, 0.0219842, 0};
Point(newp) = { 0.8407324, 0.0195838, 0};
Point(newp) = { 0.8638564, 0.0170915, 0};
Point(newp) = { 0.8868235, 0.0145051, 0};
Point(newp) = { 0.9061905, 0.0121952, 0};
Point(newp) = { 0.9225336, 0.0101138, 0};
Point(newp) = { 0.9363346, 0.0083265, 0};
Point(newp) = { 0.9479946, 0.0068038, 0};
Point(newp) = { 0.9578511, 0.0055144, 0};
Point(newp) = { 0.9661860, 0.0044240, 0};
Point(newp) = { 0.9732361, 0.0035015, 0};
Point(newp) = { 0.9792020, 0.0027211, 0};
Point(newp) = { 0.9842508, 0.0020606, 0};
Point(newp) = { 0.9885252, 0.0015014, 0};
Point(newp) = { 0.9921438, 0.0010280, 0};
Point(newp) = { 0.9952080, 0.0006271, 0};
Point(newp) = { 0.9978030, 0.0002876, 0};
Point(newp) = { 1.0000000, 0.0000000, 0};

Spline(1) = {1 : 72};
Symmetry {0, 1, 0, 0} {
    Duplicata { Curve{1}; }
}
Translate {b * Tan(a * Pi / 180), 0, b} {
    Dilate {{0, 0, 0}, {c, c, 1}} {
        Duplicata { Curve{1, 2}; }
    }
}

Curve Loop(101) = {-1, 2};
Curve Loop(102) = {-3, 4};
ThruSections(1) = {101, 102};

Sphere(2) = {0.5, 0, 0, r, 0, Pi/2, 2*Pi};
BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }

Characteristic Length{ PointsOf{ Surface{1, 2, 4}; } } = 0.01;
Mesh.MeshSizeMin = 0.01;
Mesh.MeshSizeMax = 1;
Mesh.MeshSizeFromCurvature = 1;
Mesh.MeshSizeExtendFromBoundary = 1;

Physical Surface("wing", 108) = {1, 2, 4};
Physical Surface("symmetry", 109) = {6};
Physical Surface("farfield", 110) = {5};
Physical Volume("domain", 111) = {2};
