R=2;
d=.5;

fine = 0.05;
alpha = 5;

phi = Asin(d/R);
N = 5 / fine;       // circular mesh
M = N;              // radial mesh
inner = Exp(-Log(1-d*d) / N / 2);
L = Round(-(Log(1+Exp(-2*N*Log(inner))) - Log(2)) / Log(inner));
circular = Exp(Log((Pi/2 + phi) / (Pi/2 - phi)) / L * 1.2);
radial = 1 / (3*alpha - 2);

Point(1) = {0, 0, 0};
Point(2) = {-R, 0, 0};
Point(3) = {d, Sqrt(R*R-d*d), 0};
Point(4) = {R, 0, 0};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};

Point(11) = {d, 0, 0};
Point(12) = {d+1, 0, 0};
Point(13) = {d-1, 0, 0};
Circle(11) = {12, 11, 13};

Line(21) = {2, 13};
Line(22) = {4, 12};

Line Loop(13) = {1, 2, 22, 11, -21};
Plane Surface(14) = {13};

Transfinite Line {11} = 2*N + 1 Using Progression inner;
Transfinite Line {-1} = L + 1 Using Progression circular;
Transfinite Line {-2} = 2*N-L + 1 Using Progression circular;
Transfinite Line {22, -21} = M + 1 Using Bump radial;
Transfinite Surface {14} = {2, 13, 12, 4};
Recombine Surface{14};

Extrude {0, 0, 1} {
  Surface{14};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("front") = {49};
Physical Surface("back") = {14};
Physical Surface("outer") = {32, 36};
Physical Surface("inner") = {44};
Physical Surface("symmetry") = {40, 48};
