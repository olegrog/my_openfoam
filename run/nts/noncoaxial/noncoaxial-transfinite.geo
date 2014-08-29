R=2;
d=.5;

fine = 0.05;

phi = Asin(d/R);
N = 5 / fine;       // circular mesh
M = N;              // radial mesh
circular = Exp(Log((Pi/2 + phi) / (Pi/2 - phi)) / N);
radial = Exp(Log(R) / (M - 1));

Point(1) = {0, 0, 0};
Point(2) = {R, 0, 0};
Point(3) = {-R, 0, 0};
Circle(1) = {2, 1, 3};

Point(11) = {d, 0, 0};
Point(12) = {d+1, 0, 0};
Point(13) = {d-1, 0, 0};
Circle(2) = {12, 11, 13};

Line(3) = {3, 13};
Line(4) = {2, 12};

Line Loop(13) = {3, -2, -4, 1};
Plane Surface(14) = {13};

Transfinite Line {2} = 2*N + 1;
Transfinite Line {1} = 2*N + 1 Using Progression circular;
Transfinite Line {-3, -4} = M + 1 Using Progression radial;
Transfinite Surface {14} = {3, 13, 12, 2};
Recombine Surface{14};

Extrude {0, 0, 1} {
  Surface{14};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("front") = {36};
Physical Surface("back") = {14};
Physical Surface("outer") = {35};
Physical Surface("inner") = {27};
Physical Surface("symmetry") = {23, 31};
