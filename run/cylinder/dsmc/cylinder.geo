Kn = 1;

r = .5/Kn;
R = 15*r;
Rx = R/2;
dx = (R-Rx)/2;
sw = 3.2*r;

refine_r = 8; // 6--8
refine_c = 2;
N_R = 24*2; // 24--48
N_theta = Round(N_R);

outer = .5*N_R;
inner = .8*N_R;
delta = Round(N_theta/8);

q_r = Exp(Log(refine_r*R/r)/N_R);
p_r = Exp(Log(refine_r*R/r/4)/N_R);
q_c = Exp(Log(refine_c)/N_theta);
Printf("q_r = %f, q_c = %f", q_r, q_c);
Printf("total cells = %.0f", (inner+outer)*N_theta);

Point(1) = {0, 0, 0};
Point(2) = {r, 0, 0};
Point(3) = {0, r, 0};
Point(4) = {R, 0, 0};
Point(5) = {dx, R, 0};
Point(6) = {-r, 0, 0};
Point(7) = {-Rx, 0, 0};
Point(10) = {dx, 0, 0};
Point(11) = {.5*sw, 0, 0};
Point(12) = {-sw, 0, 0};
Point(13) = {.5*sw, 1.5*sw, 0};
Point(14) = {1.5*sw, 0, 0};

Ellipse(1) = {2, 1, 3, 3};
Ellipse(2) = {3, 1, 3, 6};
Ellipse(3) = {4, 10, 5, 5};
Ellipse(4) = {5, 10, 5, 7};
Ellipse(5) = {12, 11, 13, 13};
Ellipse(6) = {13, 11, 13, 14};
Line(10) = {14, 2};
Line(11) = {14, 4};
Line(12) = {6, 12};
Line(13) = {7, 12};

Line Loop(1) = {1, 10, 6, 5, 12, 2};
Line Loop(2) = {6, 11, 3, 4, 13, 5};
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Transfinite Line {1, -2} = N_theta+1;
Transfinite Line {3, -4} = N_theta+1 Using Progression q_c;
Transfinite Line {-6} = N_theta-delta+1 Using Progression q_c;
Transfinite Line {5} = N_theta+delta+1 Using Progression q_c;
Transfinite Line {12} = inner+1 Using Bump 0.3;
Transfinite Line {-13} = outer+1 Using Progression q_r;
Transfinite Line {11} = outer+1 Using Progression Exp(Log(refine_r*R/r/10)/N_R);
Transfinite Line {-10} = inner+1 Using Progression Exp(Log(refine_r*R/r/5)/N_R);
Transfinite Surface {1} = {2, 14, 12, 6};
Transfinite Surface {2} = {14, 4, 7, 12};
Recombine Surface {1, 2};

Extrude {0, 0, 1/Kn} {
  Surface{1, 2};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1, 2};
Physical Surface("frontAndBack") = {1, 2, 45, 77};
Physical Surface("body") = {24, 28};
Physical Surface("free") = {64, 68};
Physical Surface("symmetry") = {32, 44, 60, 72};

