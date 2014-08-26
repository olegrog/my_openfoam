a0 = .3;
b0 = .7;
a1 = 1.5;
phi = 0.1 * Pi;

thick = 0.05;
fine = 0.2;

inner = 0.2 * fine;
outer = 0.05 * fine;

L = Sqrt(a0*Sin(phi)*a0*Sin(phi) + b0*Cos(phi)*b0*Cos(phi));
l = Sqrt(a0*Cos(phi)*a0*Cos(phi) + b0*Sin(phi)*b0*Sin(phi));

// inner ellipsis
Point(1) = {0, 0, 0};
Point(2) = {0, L, 0};
Point(3) = {0, -L, 0};
Point(4) = {b0*Sin(phi), b0*Cos(phi), 0};
Point(5) = {l, 0, 0};
Ellipse(1) = {2, 1, 4, 5};
Ellipse(2) = {5, 1, 4, 3};

// outer ellipsis
Point(6) = {0, 1, 0};
Point(7) = {0, -1, 0};
Point(8) = {a1, 0, 0};
Ellipse(3) = {6, 1, 8, 8};
Ellipse(4) = {8, 1, 8, 7};

Line(5) = {2, 6};
Line(6) = {3, 7};
Periodic Line {6} = {5};

Line Loop(7) = {5, 3, 4, -6, -2, -1};
Plane Surface(8) = {7};

Field[1] = Attractor;
Field[1].EdgesList = {1, 2};
Field[1].NNodesByEdge = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = outer;
Field[2].LcMax = inner;
Field[2].DistMin = thick;
Field[2].DistMax = 1;

Background Field = 2;

Extrude {0, 0, 1} {
  Surface {8};
  Layers {1};
  Recombine;
}
//Periodic Surface 19 {6, 26, 13, -30} = 31 {5, 18, -10, -17};

Physical Volume("internal") = {1};
Physical Surface("front") = {40};
Physical Surface("back") = {8};
Physical Surface("inner") = {35,39};
Physical Surface("outer") = {23,27};
Physical Surface("cyclicUpper") = {19};
Physical Surface("cyclicLower") = {31};
