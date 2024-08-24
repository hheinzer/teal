SetFactory("OpenCASCADE");
nx = 300;
ny = 240;
nz = 1;
bump = 0.03;

X = {0, 9, 14, 20, 30, 40, 55};
a = {28, 25.07355893131, 2.579601052357e1, 4.046435022819e1, 1.792461334664e1, 5.639011190988e1};
b = {0, 0.9754803562315, 8.206693007457e-1, -1.379581654948, 8.743920332081e-1, -2.010520359035};
c = {6.775070969851e-3, -1.016116352781e-1, -9.055370274339e-2, 1.945884504128e-2, -5.567361123058e-2, 1.644919857549e-2};
d = {-2.124527775800e-3, 1.889794677828e-3, 1.626510569859e-3, -2.070318932190e-4, 6.277731764683e-4, 2.674976141766e-5};

For i In {0: 5}
    For x In {X[i]: X[i + 1] - 1}
        Point(newp) = {x / 28, Max(0, Min(a[i] + b[i] * x + c[i] * x^2 + d[i] * x^3, 28)) / 28, 0};
    EndFor
EndFor
For i In {0: 5}
    For x In {X[i]: X[i + 1] - 1}
        Point(newp) = {(252 - x) / 28, Max(0, Min(a[i] + b[i] * x + c[i] * x^2 + d[i] * x^3, 28)) / 28, 0};
    EndFor
EndFor
Point(newp) = {9, 3.035, 0};
Point(newp) = {(252 - 54) / 28, 3.035, 0};
Point(newp) = {54 / 28, 3.035, 0};
Point(newp) = {0, 3.035, 0};

Spline(1) = {1: 55};
Line(2) = {55, 110};
Spline(3) = {110: 56};
Line(4) = {56, 111};
Line(5) = {111, 112};
Line(6) = {112, 110};
Line(7) = {112, 113};
Line(8) = {113, 55};
Line(9) = {113, 114};
Line(10) = {114, 1};

Curve Loop(1) = {1, -8, 9, 10};
Plane Surface(1) = {1};
Curve Loop(2) = {2, -6, 7, 8};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 4, 5, 6};
Plane Surface(3) = {3};

Extrude {0, 0, 4.5} {
  Surface{1}; Surface{2}; Surface{3}; Layers {nz}; Recombine;
}

Physical Surface("inlet", 29) = {7};
Physical Surface("outlet", 30) = {14};
Physical Surface("wall", 31) = {4, 9, 13, 15, 11, 6};
Physical Surface("back", 32) = {1, 2, 3};
Physical Surface("front", 33) = {8, 12, 16};
Physical Volume("domain", 34) = {1, 2, 3};

Transfinite Curve {10, 8, 6, 4} = (ny + 1) Using Bump bump;
Transfinite Curve {1, 3, 5, 9} = ((54 / 28) * (nx + 1) / 9 + 1) Using Progression 1;
Transfinite Curve {2, 7} = ((9 - 2 * 54 / 28) * (nx + 1) / 9 + 1) Using Progression 1;
Transfinite Surface {1, 2, 3};
Recombine Surface {1, 2, 3};
