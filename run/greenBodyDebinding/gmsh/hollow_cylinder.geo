SetFactory("OpenCASCADE");

R = 5e-3;
r = 0.8*R;
N = 15; // ratio to cell size

Cylinder(1) = {0, 0, 0, 0, 0, 2*R, R};
Cylinder(2) = {0, 0, 0, 0, 0, 2*R, r};
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

MeshSize{:} = R/N;

Physical Volume("internal") = {3};
Physical Surface("all") = Boundary{ Volume{3}; };

