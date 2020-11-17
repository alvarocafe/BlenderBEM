// Gmsh project created on Sun Aug 28 10:42:10 2016
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {1, 0, 0, 1.0};
Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Line("left") = {1};
Physical Line("down") = {2};
Physical Line("right") = {3};
Physical Line("upper") = {4};
Physical Surface("surface") = {6};
