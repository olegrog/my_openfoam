a0=.3;
b0=.7;
a1=1.5;

thick = 0;  //0
fluid = 0.03; //0.04
bound = 0.0002;//0.005

Point(1) = {0,  0, 0};
Point(2) = {a0,  0, 0};
Point(3) = {0 , b0, 0};
Point(4) = {a1,  0, 0};
Point(5) = {0 ,  1, 0};
Ellipse(1) = {2,1,3,3};
Ellipse(2) = {4,1,5,5};
Line(3)   = {2,4};
Line(4)   = {3,5};

Line Loop(5) = {2, -4, -1, 3};
Plane Surface(6) = {5};

Field[1] = Attractor;
Point(10) = {0, b0+0.02, 0};
Line(10)  = {3, 10};
Field[1].EdgesList = {1,10};
Field[1].NNodesByEdge = 1000;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = bound;
Field[2].LcMax = fluid;
Field[2].DistMin = thick;
Field[2].DistMax = .05; //2.5

Background Field = 2;


Extrude {0, 0, 1} {
  Surface{6};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("front") = {32};
Physical Surface("back") = {6};
Physical Surface("inner") = {27};
Physical Surface("outer") = {19};
Physical Surface("symmetry") = {23, 31};
