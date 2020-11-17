SetFactory("OpenCASCADE");

R = 5e-3;
N = 10; // ratio to cell size

Sphere(1) = {0, 0, 0, R};

MeshSize{:} = R/N;

Physical Volume("internal") = {1};
Physical Surface("all") = Boundary{ Volume{1}; };

