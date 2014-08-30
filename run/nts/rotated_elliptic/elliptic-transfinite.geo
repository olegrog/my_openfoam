a0 = .3;
b0 = .7;
a1 = 1.5;
// 0 < betta < Pi/2
betta = 0.1 * Pi;

fine = 0.05;

/***** inner ellipsis *****/
t = Atan2(b0*Sin(betta), a0*Cos(betta));
L = Sqrt(b0*Cos(t)*b0*Cos(t) + a0*Sin(t)*a0*Sin(t));
Point(1) = {0, 0, 0};
Point(2) = {0, L, 0};
Point(3) = {b0*Sin(betta), b0*Cos(betta), 0};
Point(4) = {a0*Cos(betta), -a0*Sin(betta), 0};
Point(5) = {0, -L, 0};
Ellipse(1) = {2, 1, 3, 3};
Ellipse(2) = {3, 1, 3, 4};
Ellipse(3) = {4, 1, 3, 5};

/***** outer ellipsis *****/
u = Pi - Atan2(a1*Cos(betta), -Sin(betta));
v = Atan2(a1*Sin(betta), Cos(betta));
Point(6) = {0, 1, 0};
Point(7) = {a1*Cos(u), Sin(u), 0};
Point(8) = {a1, 0, 0};
Point(9) = {a1*Cos(v), -Sin(v), 0};
Point(10) = {0, -1, 0};
Ellipse(4) = {6, 1, 8, 7};
Ellipse(5) = {7, 1, 8, 8};
Ellipse(6) = {8, 1, 8, 9};
Ellipse(7) = {9, 1, 8, 10};

/***** other lines *****/
Line(10) = {2, 6};
Line(11) = {5, 10};
Periodic Line {11} = {10};

Line Loop(7) = {10, 4, 5, 6, 7, -11, -3, -2, -1};
Plane Surface(8) = {7};

N = 5 / fine;           // circular mesh
M = Round(2*t*N/Pi);
L = N;                  // radial mesh
outer = Exp(Log(a1*(M*(Pi/2-v)) / ((N-M)*(Pi/2-u))) / N);
N1 = Round(Log(1-u/(u+v)*(1-Exp(N*Log(outer))))/Log(outer));
inner = Exp(Log(b0/a0)/N);
Printf("%g %g %g", N, N1, outer);
radial = Exp(Log(1/a0) / (L - 1));

Transfinite Line {1} = M + 1 Using Progression 1/inner;
Transfinite Line {2} = N + 1 Using Progression inner;
Transfinite Line {3} = N - M + 1 Using Progression 1/inner;
Transfinite Line {4} = M + 1;
Transfinite Line {5} = N1 + 1 Using Progression outer;
Transfinite Line {6} = N - N1 + 1 Using Progression outer;
Transfinite Line {7} = N - M + 1;
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
Physical Surface("front") = {58};
Physical Surface("inner") = {49, 53, 57};
Physical Surface("outer") = {29, 33, 37, 41};
Physical Surface("cyclicUpper") = {25};
Physical Surface("cyclicLower") = {45};
