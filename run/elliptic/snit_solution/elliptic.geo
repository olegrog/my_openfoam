Kn = 0.02;
a0 = .3;
b0 = .7;
a1 = 1.5;
b1 = 1;

refine_i = 10; //40
refine_o = refine_i / 2; // /4

L = 40;
N = 1.2*L;

d_i = 0.2;
d_o = 0.1; //0.05
c_i = 2.5; //4
c_o = 2; //3

inner = Exp(Log(2*b0/a0)/N);
outer = Exp(Log(2*b0/a0)/N);
i_axes = Exp(Log(refine_i)/(d_i*L*c_i));
o_axes = Exp(Log(refine_o)/(d_o*L*c_o));

Point(1) = {0,   0, 0};
Point(2) = {a0,  0, 0};
Point(3) = {0 , b0, 0};
Point(4) = {a1,  0, 0};
Point(5) = {0 , b1, 0};

Point(6) = {a0 + d_i*(a1-a0), 0, 0};
Point(7) = {a1 - d_o*(a1-a0), 0, 0};
Point(8) = {0, b0 + d_i*(b1-b0), 0};
Point(9) = {0, b1 - d_o*(b1-b0), 0};

Ellipse(1) = {2, 1, 3, 3};
Ellipse(2) = {4, 1, 5, 5};
Line(3) = {2, 6};
Line(4) = {6, 7};
Line(5) = {7, 4};
Line(6) = {3, 8};
Line(7) = {8, 9};
Line(8) = {9, 5};

Line Loop(5) = {2, -8, -7, -6, -1, 3, 4, 5};
Plane Surface(6) = {5};

Printf("N = %f, inner = %f", N, inner);
Printf("N_uniform = %f", (L+1)*(1-d_o-d_i));
Printf("q^n_inner = %f, q^n_outer = %f", Exp(Log(i_axes)*c_i*(L+1)*d_i), Exp(Log(o_axes)*c_o*(L+1)*d_o));
Transfinite Line {1} = N+1 Using Progression 1/inner;
Transfinite Line {2} = N+1 Using Progression 1/outer;
Transfinite Line {4, 7} = (L+1)*(1-d_o-d_i);
Transfinite Line {3, 6} = c_i*(L+1)*d_i Using Progression i_axes;
Transfinite Line {-5, -8} = c_o*(L+1)*d_o Using Progression o_axes;
Transfinite Surface {6} = {2, 4, 5, 3};
Recombine Surface {6};

Extrude {0, 0, 1} {
  Surface{6};
  Layers{1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("front") = {50};
Physical Surface("back") = {6};
Physical Surface("inner") = {37};
Physical Surface("outer") = {21};
Physical Surface("left") = {25,29,33};
Physical Surface("bottom") = {41,45,49};
