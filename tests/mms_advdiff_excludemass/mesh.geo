res=${res}; // resolution
Point(1) = {0.1,-0.3,0,res};
Point(2) = {0.6,-0.3,0,res};
Point(3) = {0.6,0.1,0,res};
Point(4) = {0.1,0.1,0,res};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
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
