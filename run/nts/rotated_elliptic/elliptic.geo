a0 = .3;
b0 = .7;
a1 = 1.5;
// 0 <= phi < Pi
phi = 0.1 * Pi;

thick = 0.05;
fine = 0.2;

inner = 0.2 * fine;
outer = 0.05 * fine;

/***** inner ellipsis *****/
t = Atan2(a0*Sin(phi), b0*Cos(phi));
l = Sqrt(a0*Cos(t)*a0*Cos(t) + b0*Sin(t)*b0*Sin(t));
t = Atan2(b0*Sin(phi), a0*Cos(phi));
L = Sqrt(b0*Cos(t)*b0*Cos(t) + a0*Sin(t)*a0*Sin(t));
Point(1) = {0, 0, 0};
Point(2) = {0, L, 0};
Point(3) = {b0*Sin(phi), b0*Cos(phi), 0};
Point(4) = {l, 0, 0};
Point(5) = {0, -L, 0};
If (phi != 0)
    Ellipse(1) = {2, 1, 3, 3};
EndIf
If (phi != 0.5 * Pi)
    Ellipse(2) = {3, 1, 3, 4};
EndIf
If (phi != Pi)
    Ellipse(3) = {4, 1, 3, 5};
EndIf

/***** outer ellipsis *****/
Point(6) = {0, 1, 0};
Point(7) = {a1, 0, 0};
Point(8) = {0, -1, 0};
Ellipse(4) = {6, 1, 7, 7};
Ellipse(5) = {7, 1, 7, 8};

/***** other lines *****/
Line(6) = {2, 6};
Line(7) = {5, 8};
Periodic Line {7} = {6};

Line Loop(7) = {6, 4, 5, -7, -3, -2, -1};
Plane Surface(8) = {7};

Field[1] = Attractor;
Field[1].EdgesList = {1, 2};
Field[1].NNodesByEdge = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = outer;
Field[2].LcMax = inner;
Field[2].DistMin = thick;
Field[2].DistMax = 1;

Background Field = 2;

Extrude {0, 0, 1} {
  Surface {8};
  Layers {1};
  Recombine;
}

Physical Volume("internal") = {1};
Physical Surface("back") = {8};
If (phi == 0 || phi == 0.5 * Pi)
    Physical Surface("front") = {40};
    Physical Surface("inner") = {35, 39};
    Physical Surface("outer") = {23, 27};
    Physical Surface("cyclicUpper") = {19};
    Physical Surface("cyclicLower") = {31};
    Abort;
EndIf
Physical Surface("front") = {45};
Physical Surface("inner") = {36, 40, 44};
Physical Surface("outer") = {24, 28};
Physical Surface("cyclicUpper") = {20};
Physical Surface("cyclicLower") = {32};
