ncells = ${ncells};
Point(1) = {-2, -2, 0, 1.0};
Extrude {1, 0, 0} {
  Point{1};
}
Extrude {1, 0, 0} {
  Point{2};
}
Extrude {1, 0, 0} {
  Point{3};
}
Extrude {1, 0, 0} {
  Point{4};
}
Extrude {0, 1, 0} {
  Line{1, 2, 3, 4};
}
Extrude {0, 1, 0} {
  Line{5, 13, 9, 17};
}
Extrude {0, 1, 0} {
  Line{21, 29, 25, 33};
}
Extrude {0, 1, 0} {
  Line{37, 41, 45, 49};
}
Transfinite Line {1}  = ncells+1 Using Progression 1;
Transfinite Line {2}  = ncells+1 Using Progression 1;
Transfinite Line {3}  = ncells+1 Using Progression 1;
Transfinite Line {4}  = ncells+1 Using Progression 1;
Transfinite Line {5}  = ncells+1 Using Progression 1;
Transfinite Line {9}  = ncells+1 Using Progression 1;
Transfinite Line {13} = ncells+1 Using Progression 1;
Transfinite Line {17} = ncells+1 Using Progression 1;
Transfinite Line {21} = ncells+1 Using Progression 1;
Transfinite Line {25} = ncells+1 Using Progression 1;
Transfinite Line {29} = ncells+1 Using Progression 1;
Transfinite Line {33} = ncells+1 Using Progression 1;
Transfinite Line {37} = ncells+1 Using Progression 1;
Transfinite Line {41} = ncells+1 Using Progression 1;
Transfinite Line {45} = ncells+1 Using Progression 1;
Transfinite Line {49} = ncells+1 Using Progression 1;
Transfinite Line {53} = ncells+1 Using Progression 1;
Transfinite Line {57} = ncells+1 Using Progression 1;
Transfinite Line {61} = ncells+1 Using Progression 1;
Transfinite Line {65} = ncells+1 Using Progression 1;
Transfinite Line {6}  = ncells+1 Using Progression 1;
Transfinite Line {22} = ncells+1 Using Progression 1;
Transfinite Line {38} = ncells+1 Using Progression 1;
Transfinite Line {54} = ncells+1 Using Progression 1;
Transfinite Line {7}  = ncells+1 Using Progression 1;
Transfinite Line {23} = ncells+1 Using Progression 1;
Transfinite Line {39} = ncells+1 Using Progression 1;
Transfinite Line {55} = ncells+1 Using Progression 1;
Transfinite Line {11} = ncells+1 Using Progression 1;
Transfinite Line {26} = ncells+1 Using Progression 1;
Transfinite Line {43} = ncells+1 Using Progression 1;
Transfinite Line {59} = ncells+1 Using Progression 1;
Transfinite Line {15} = ncells+1 Using Progression 1;
Transfinite Line {27} = ncells+1 Using Progression 1;
Transfinite Line {47} = ncells+1 Using Progression 1;
Transfinite Line {63} = ncells+1 Using Progression 1;
Transfinite Line {19} = ncells+1 Using Progression 1;
Transfinite Line {35} = ncells+1 Using Progression 1;
Transfinite Line {51} = ncells+1 Using Progression 1;
Transfinite Line {67} = ncells+1 Using Progression 1;

Transfinite Surface {8}  Alternated;
Transfinite Surface {12} Alternated;
Transfinite Surface {16} Alternated;
Transfinite Surface {20} Alternated;
Transfinite Surface {36} Alternated;
Transfinite Surface {28} Alternated;
Transfinite Surface {32} Alternated;
Transfinite Surface {24} Alternated;
Transfinite Surface {40} Alternated;
Transfinite Surface {44} Alternated;
Transfinite Surface {48} Alternated;
Transfinite Surface {52} Alternated;
Transfinite Surface {68} Alternated;
Transfinite Surface {64} Alternated;
Transfinite Surface {60} Alternated;
Transfinite Surface {56} Alternated;


Physical Line(1) = {23, 39};
Physical Line(2) = {27, 47};
Physical Line(3) = {9, 13};
Physical Line(4) = {41, 45};

// Center
Physical Surface(1) = {32};
Physical Surface(2) = {28};
Physical Surface(3) = {44};
Physical Surface(4) = {48};

// Left
Physical Surface(10) = {24};
Physical Surface(15) = {40};

// Right
Physical Surface(20) = {36};
Physical Surface(25) = {52};

// Bottom
Physical Surface(30) = {12};
Physical Surface(35) = {16};

// Top
Physical Surface(40) = {60};
Physical Surface(45) = {64};

// Bottom Left
Physical Surface(50) = {8};

// Bottom Right
Physical Surface(60) = {20};

// Top Left
Physical Surface(70) = {56};

// Top Right
Physical Surface(80) = {68};
