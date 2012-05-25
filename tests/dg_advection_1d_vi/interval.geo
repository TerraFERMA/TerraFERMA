Point(1) = {0, 0, 0, 0.1};
Extrude {0.5, 0, 0} {
  Point{1}; Layers{5};
}
Extrude {4.5, 0, 0} {
  Point{2}; Layers{45};
}
Physical Point(1) = {1};
Physical Point(2) = {3};
Physical Line(1) = {1};
Physical Line(2) = {2};
