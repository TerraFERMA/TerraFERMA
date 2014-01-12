nhpoints=${nhpoints}; // resolution x
nvpoints=${nvpoints}; // resolution y
Point(1) = {0.1,-0.3,0};
Point(2) = {0.6,-0.3,0};
Point(3) = {0.6,0.1,0};
Point(4) = {0.1,0.1,0};
Line(1) = {1,2};
Transfinite Line{1} = nhpoints;
Line(2) = {2,3};
Transfinite Line{2} = nvpoints;
Line(3) = {3,4};
Transfinite Line{3} = nhpoints;
Line(4) = {4,1};
Transfinite Line{4} = nvpoints;
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Transfinite Surface{6} = {1,2,3,4} Alternate;
// y const
Physical Line(7) = {1};
// x const
Physical Line(8) = {2};
// y const
Physical Line(9) = {3};
// x const
Physical Line(10) = {4};
// z const
Physical Surface(11) = {6};
