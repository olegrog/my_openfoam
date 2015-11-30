Kn = 0.01;
D = 1;
L1 = 0.15*D;
L2 = D;
d = 1;
theta = Pi/50;

Nx = 20;

h_max = (D/2) / Nx;
R = 0.9; Nx = 45;
h_min = 0.1*Kn;
h_min = h_min < h_max/4 ? h_min : h_max/4;
l = (1+Fabs(Log(Kn))/3)*Kn;
l = l < D/6 ? l : D/6;
l = l > 6*h_max ? l : 6*h_max;

Ny_fl = Round((D/2-l)/(D/2)*Nx);
r = 0.3;
p1 = Exp(Log(L1/D*r/2)/Nx);
p2 = Exp(Log(L2/D*r/2)/Nx);

Printf("h_max = %f, h_min = %f", h_max, h_min);
a = h_min/h_max;
b = l/h_max;
q = h_max > h_min ? (1-b)/(a-b) : 0;
Printf("q = %f", q);
Ny_kn = 1 + Round(Log(a)/Log(q));
q = Exp(Log(a)/(Ny_kn-1));

Printf("l = %f, D/2 = %f", l, D/2);
Printf("Size: %f x %f", Nx*(1.5*D+L1+L2)/D*2, Ny_fl+Ny_kn);
Printf("q = %f, p1 = %f, p2 = %f", q, p1, p2);
If (q < 0.5)
    //Error("q is too small!");
EndIf

Point(1) =  { 0,           0, 0};
Point(2) =  { D,           0, 0};
Point(3) =  { D+L1,        0, 0};
Point(4) =  { 2*D+L1,      0, 0};
Point(5) =  { 2*D+L1+L2,   0, 0};
Point(6) =  { 3*D+L1+L2,   0, 0};

Point(11) =  { 0,           l, 0};
Point(12) =  { D,           l, 0};
Point(13) =  { D+L1,        l, 0};
Point(14) =  { 2*D+L1,      l, 0};
Point(15) =  { 2*D+L1+L2,   l, 0};
Point(16) =  { 3*D+L1+L2,   l, 0};

Point(21) =  { 0,           D/2, 0};
Point(22) =  { D,           D/2, 0};
Point(23) =  { D+L1,        D/2, 0};
Point(24) =  { 2*D+L1,      D/2, 0};
Point(25) =  { 2*D+L1+L2,   D/2, 0};
Point(26) =  { 3*D+L1+L2,   D/2, 0};

Line(1) = {1, 11};
Line(2) = {2, 12};
Line(3) = {3, 13};
Line(4) = {4, 14};
Line(5) = {5, 15};
Line(6) = {6, 16};

Line(11) = {11, 21};
Line(12) = {12, 22};
Line(13) = {13, 23};
Line(14) = {14, 24};
Line(15) = {15, 25};
Line(16) = {16, 26};

Line(21) = {1, 2};
Line(22) = {2, 3};
Line(23) = {3, 4};
Line(24) = {4, 5};
Line(25) = {5, 6};

Line(31) = {11, 12};
Line(32) = {12, 13};
Line(33) = {13, 14};
Line(34) = {14, 15};
Line(35) = {15, 16};

Line(41) = {21, 22};
Line(42) = {22, 23};
Line(43) = {23, 24};
Line(44) = {24, 25};
Line(45) = {25, 26};

Line Loop(1) = {1, 11, 41, -12, -2, -21};
Line Loop(2) = {2, 12, 42, -13, -3, -22};
Line Loop(3) = {3, 13, 43, -14, -4, -23};
Line Loop(4) = {4, 14, 44, -15, -5, -24};
Line Loop(5) = {5, 15, 45, -16, -6, -25};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

Transfinite Line{11} = Ny_fl + 1 Using Progression R;
Transfinite Line{12} = Ny_fl + 1 Using Progression R;
Transfinite Line{13} = Ny_fl + 1 Using Progression R;
Transfinite Line{14} = Ny_fl + 1 Using Progression R;
Transfinite Line{15} = Ny_fl + 1 Using Progression R;
Transfinite Line{16} = Ny_fl + 1 Using Progression R;

Transfinite Line{-1} = Ny_kn + 1 Using Progression q;
Transfinite Line{-2} = Ny_kn + 1 Using Progression q;
Transfinite Line{-3} = Ny_kn + 1 Using Progression q;
Transfinite Line{-4} = Ny_kn + 1 Using Progression q;
Transfinite Line{-5} = Ny_kn + 1 Using Progression q;
Transfinite Line{-6} = Ny_kn + 1 Using Progression q;

Transfinite Line{21} = Nx + 1 Using Progression p1;
Transfinite Line{22} = Nx + 1 Using Bump r;
Transfinite Line{23} = Nx + 1 Using Bump r*L1/D;
Transfinite Line{24} = Nx + 1 Using Bump r;
Transfinite Line{25} = Nx + 1 Using Progression 1/p2;

Transfinite Line{31} = Nx + 1 Using Progression p1;
Transfinite Line{32} = Nx + 1 Using Bump r;
Transfinite Line{33} = Nx + 1 Using Bump r*L1/D;
Transfinite Line{34} = Nx + 1 Using Bump r;
Transfinite Line{35} = Nx + 1 Using Progression 1/p2;

Transfinite Line{41} = Nx + 1 Using Progression p1;
Transfinite Line{42} = Nx + 1 Using Bump r;
Transfinite Line{43} = Nx + 1 Using Bump r*L1/D;
Transfinite Line{44} = Nx + 1 Using Bump r;
Transfinite Line{45} = Nx + 1 Using Progression 1/p2;

Transfinite Surface{1} = {1, 21, 22, 2};
Transfinite Surface{2} = {2, 22, 23, 3};
Transfinite Surface{3} = {3, 23, 24, 4};
Transfinite Surface{4} = {4, 24, 25, 5};
Transfinite Surface{5} = {5, 25, 26, 6};
Recombine Surface{1, 2, 3, 4, 5};

Rotate {{1, 0, 0}, {0, D/2, 0}, -theta/2} {
    Surface{1, 2, 3, 4, 5};
}
Extrude {{1, 0, 0}, {0, D/2, 0}, theta} {
    Surface{1, 2, 3, 4, 5};
    Layers{1};
    Recombine;
}

Physical Surface("hot") = {71, 179};
Physical Surface("cold") = {125};
Physical Surface("short") = {98};
Physical Surface("long") = {152};
Physical Surface("left") = {56, 59};
Physical Surface("right") = {171, 175};
Physical Surface("wedge1") = {1, 2, 3, 4, 5};
Physical Surface("wedge2") = {72, 99, 126, 153, 180};

Physical Volume("coldV") = {3};
Physical Volume("hotV") = {1, 5};
Physical Volume("shortV") = {2};
Physical Volume("longV") = {4};
