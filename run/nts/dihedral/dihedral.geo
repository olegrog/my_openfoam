Point(1) = {0, 0, 0};
Point(2) = {5, 0, 0};
Point(3) = {4, 3, 0};
Point(4) = {1, 0, 0};
Point(5) = {.8, .6, 0};

Line(1) = {2, 4};
Line(2) = {3, 5};
Circle(3) = {2, 1, 3};
Circle(4) = {4, 1, 5};

Line Loop(5) = {2, -4, -1, 3};
Plane Surface(6) = {5};

Field[1] = Attractor;
Field[1].NodesList = {4};
Field[1].NNodesByEdge = 50;

k = 4;
p = .02;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = p;
Field[2].LcMax = p*k;
Field[2].DistMin = 0;
Field[2].DistMax = k;

Background Field = 2;

Extrude {0, 0, 1} {
  Surface{6};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("front") = {28};
Physical Surface("back") = {6};
Physical Surface("bottom") = {23};
Physical Surface("top") = {15};
Physical Surface("inlet") = {19};
Physical Surface("outlet") = {27};
