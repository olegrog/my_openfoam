r=2;
d=.5;

thick = 0.05;
fine = 0.2;

inner = 0.2 * fine;
outer = 0.05 * fine;

Point(1) = {0,0,0};
Point(2) = {r , 0, 0};
Point(3) = {-r, 0, 0};
Point(4) = {0 , r, 0};
Circle(1) = {2,1,4};
Circle(2) = {4,1,3};

Point(11) = {d,0,0};
Point(12) = {d+1 , 0, 0};
Point(13) = {d-1, 0, 0};
Point(14) = {d , 1, 0};
Circle(11) = {12,11,14};
Circle(12) = {14,11,13};

Line(3)    = {3,13};
Line(4)    = {2,12};

Line Loop(13) = {2, 3, -12, -11, -4, 1};
Plane Surface(14) = {13};

Field[1] = Attractor;
Field[1].EdgesList = {11,12};
Field[1].NNodesByEdge = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = outer;
Field[2].LcMax = inner;
Field[2].DistMin = thick;
Field[2].DistMax = 1;

Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;

Extrude {0, 0, 1} {
  Surface{14};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("front") = {46};
Physical Surface("back") = {14};
Physical Surface("outer") = {25, 45};
Physical Surface("inner") = {37, 33};
Physical Surface("symmetry") = {29, 41};
Coherence;
