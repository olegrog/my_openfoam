Point(1) = {-1, 0, 0};
Point(2) = {1 , 0, 0};
Point(3) = {0 , 2, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};

Line Loop(1) = {1, 2, 3};
Plane Surface(2) = {1};

//Recombine Surface {2};

Extrude {0, 0, 1} {
  Surface{2};
  Layers{1};
  Recombine;
}

lc = 0.025;
R = 100;

Field[1] = Distance;
Field[1].PointsList = {1,2,3};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc / R;
Field[2].SizeMax = lc;
Field[2].DistMin = 0;
Field[2].DistMax = 0.5;
Field[2].Sigmoid = 0;
Background Field = 2;

Physical Volume("internal") = {1};
Physical Surface("frontAndBack") = {2, 20};
Physical Surface("side") = {15, 19};
Physical Surface("bottom") = {11};
