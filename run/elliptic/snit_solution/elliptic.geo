a0 = .3;
b0 = .7;
a1 = 1.5;

fine = 0.1;
refine = 0.03;

N = 5 / fine;
L = N;
inner = Exp(Log(b0/a0)/N);

Point(1) = {0,   0, 0};
Point(2) = {a0,  0, 0};
Point(3) = {0 , b0, 0};
Point(4) = {a1,  0, 0};
Point(5) = {0 ,  1, 0};
Ellipse(1) = {2, 1, 3, 3};
Ellipse(2) = {4, 1, 5, 5};
Line(3) = {2, 4};
Line(4) = {3, 5};

Line Loop(5) = {2, -4, -1, 3};
Plane Surface(6) = {5};

Printf("N = %f, inner = %f", N, inner);
Transfinite Line {1} = N + 1 Using Progression 1/inner;
Transfinite Line {2} = N + 1;
Transfinite Line {3, 4} = L + 1 Using Bump refine;
Transfinite Surface {6} = {2, 4, 5, 3};
Recombine Surface {6};

Extrude {0, 0, 1} {
  Surface{6};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("front") = {28};
Physical Surface("back") = {6};
Physical Surface("inner") = {23};
Physical Surface("outer") = {15};
Physical Surface("left") = {19};
Physical Surface("bottom") = {27};
