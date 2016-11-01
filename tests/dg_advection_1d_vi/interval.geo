Point(1) = {0, 0, 0, 0.1};
Extrude {0.5, 0, 0} {
  Point{1}; Layers{5};
}
Extrude {4.0, 0, 0} {
  Point{2}; Layers{40};
}
Extrude {0.5, 0, 0} {
  Point{3}; Layers{5};
}
Physical Point(1) = {1};
Physical Point(2) = {4};
Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
