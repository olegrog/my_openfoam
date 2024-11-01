a0 = .3;
b0 = .7;
a1 = 1.5;
// beta = 0

fine = 0.05;

/***** inner ellipsis *****/
Point(1) = {0, 0, 0};
Point(2) = {0, b0, 0};
Point(3) = {a0, 0, 0};
Point(5) = {0, -b0, 0};
Ellipse(1) = {2, 1, 3, 3};
Ellipse(3) = {3, 1, 3, 5};

/***** outer ellipsis *****/
Point(6) = {0, 1, 0};
Point(8) = {a1, 0, 0};
Point(10) = {0, -1, 0};
Ellipse(4) = {6, 1, 8, 8};
Ellipse(7) = {8, 1, 8, 10};

/***** other lines *****/
Line(10) = {2, 6};
Line(11) = {5, 10};
Periodic Line {11} = {10};

Line Loop(7) = {10, 4, 7, -11, -3, -1};
Plane Surface(8) = {7};

N = 5 / fine;           // circular mesh
L = N;                  // radial mesh
inner = Exp(Log(b0/a0)/N);
radial = Exp(Log(1/a0) / (L - 1));

Transfinite Line {1} = N + 1 Using Progression inner;
Transfinite Line {3} = N + 1 Using Progression 1/inner;
Transfinite Line {4} = N + 1;
Transfinite Line {7} = N + 1;
Transfinite Line {10, 11} = L + 1 Using Progression radial;
Transfinite Surface {8} = {2, 6, 10, 5};
Recombine Surface{8};

Extrude {0, 0, 1} {
  Surface {8};
  Layers {1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("back") = {8};
Physical Surface("front") = {43};
Physical Surface("inner") = {38, 42};
Physical Surface("outer") = {26, 30};
Physical Surface("cyclicUpper") = {22};
Physical Surface("cyclicLower") = {34};
