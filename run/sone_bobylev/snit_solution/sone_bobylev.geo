Kn = 1;
Nx = 30;
L = 1;

threshold = 0.1;
If (Kn > threshold || Kn == 0)
    Kn = threshold;
EndIf
Kn = 0.1;

h_max = (L/2) / Nx;
h_min = 0.02*Kn;
H = Exp(2.5*Log(-Log(Kn)/Log(10)))*Kn;
l = L/2 > H ? H : L/2;

Ny_fl = L/2 > H ? Round((L/2-l)/(L/2)*Nx) : Nx;

a = h_min/h_max;
b = l/h_max;
q = h_max > h_min ? (1-b)/(a-b) : 0;
Ny_kn = 1 + Round(Log(a)/Log(q));
q = Exp(Log(a)/(Ny_kn-1));

Printf("l = %f, L/2 = %f", l, L/2);
Printf("Size: %f x %f", Nx, Ny_fl+Ny_kn);
Printf("q = %f", q);
If (q < 0.5)
    Error("q is too small!");
EndIf

Point(1) =  { 0,   0,   0};
Point(2) =  { 0,   L/2, 0};
Point(12) = { 0,   l,   0};
Point(3) =  { L/2, L/2, 0};
Point(13) = { L/2, l,   0};
Point(4) =  { L/2, 0,   0};

Line(1) = {1, 12};
Line(11) = {12, 2};
Line(2) = {2, 3};
Line(13) = {3, 13};
Line(3) = {13, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 11, 2, 13, 3, 4};
Plane Surface(1) = {1};

Transfinite Line{-11} = Ny_fl + 1;
Transfinite Line{-1} = Ny_kn + 1 Using Progression q;
Transfinite Line{2} = Nx + 1;
Transfinite Line{3} = Ny_kn + 1 Using Progression q;
Transfinite Line{13} = Ny_fl + 1;
Transfinite Line{4} = Nx + 1;
Transfinite Surface{1} = {1, 2, 3, 4};
Recombine Surface{1};

Extrude {0, 0, L/10} {
	Surface{1};
	Layers{1};
	Recombine;
}

Physical Surface("bottom") = {44};
Physical Surface("top") = {32};
Physical Surface("left") = {24, 28};
Physical Surface("right") = {36, 40};
Physical Surface("frontAndBack") = {1, 45};
Physical Volume("volume") = {1};
