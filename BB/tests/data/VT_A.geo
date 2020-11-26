
// Defining the points
Point(1) = {-7.569397566060481, 0, 0, 10.0};
Point(2) = {0, 0, 0,10.0};
Point(3) = {7.569397566060481,0,0,10.0};
// Defining the circles
Circle(1) = {1, 2, 3};
Circle(2) = {3, 2, 1};
// Definig the surface area of the glotis
Line Loop(1) = {-2, -1};
Plane Surface(1) = {1};

// Defining the points
Point(4) = {-7.569397566060481, 0, 10, 10.0};
Point(5) = {0, 0, 10,10.0};
Point(6) = {7.569397566060481,0,10,10.0};
// Defining the circles
Circle(5) = {4, 5, 6};
Circle(6) = {6, 5, 4};
// Defining the surface
Line Loop(2) = {6, 5};


// Defining the lines between circle 5 and 4
Line(7) = {1, 4};
Line(8) = {3, 6};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(3) = {1, 8, -5, -7};
Ruled Surface(2) = {3};
Line Loop(4) = {2, 7, -6, -8};
Ruled Surface(3) = {4};
			
// Defining the points
Point(7) = {-9.440697438826296, 0, 20, 10.0};
Point(8) = {0, 0, 20,10.0};
Point(9) = {9.440697438826296,0,20,10.0};
// Defining the circles
Circle(9) = {7, 8, 9};
Circle(10) = {9, 8, 7};
// Defining the surface
Line Loop(5) = {10, 9};


// Defining the lines between circle 9 and 8
Line(11) = {4, 7};
Line(12) = {6, 9};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(6) = {5, 12, -9, -11};
Ruled Surface(4) = {6};
Line Loop(7) = {6, 11, -10, -12};
Ruled Surface(5) = {7};
			
// Defining the points
Point(10) = {-6.90988298942671, 0, 30, 10.0};
Point(11) = {0, 0, 30,10.0};
Point(12) = {6.90988298942671,0,30,10.0};
// Defining the circles
Circle(13) = {10, 11, 12};
Circle(14) = {12, 11, 10};
// Defining the surface
Line Loop(8) = {14, 13};


// Defining the lines between circle 13 and 12
Line(15) = {7, 10};
Line(16) = {9, 12};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(9) = {9, 16, -13, -15};
Ruled Surface(6) = {9};
Line Loop(10) = {10, 15, -14, -16};
Ruled Surface(7) = {10};
			
// Defining the points
Point(13) = {-5.046265044040321, 0, 40, 10.0};
Point(14) = {0, 0, 40,10.0};
Point(15) = {5.046265044040321,0,40,10.0};
// Defining the circles
Circle(17) = {13, 14, 15};
Circle(18) = {15, 14, 13};
// Defining the surface
Line Loop(11) = {18, 17};


// Defining the lines between circle 17 and 16
Line(19) = {10, 13};
Line(20) = {12, 15};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(12) = {13, 20, -17, -19};
Ruled Surface(8) = {12};
Line Loop(13) = {14, 19, -18, -20};
Ruled Surface(9) = {13};
			
// Defining the points
Point(16) = {-6.432750982580687, 0, 50, 10.0};
Point(17) = {0, 0, 50,10.0};
Point(18) = {6.432750982580687,0,50,10.0};
// Defining the circles
Circle(21) = {16, 17, 18};
Circle(22) = {18, 17, 16};
// Defining the surface
Line Loop(14) = {22, 21};


// Defining the lines between circle 21 and 20
Line(23) = {13, 16};
Line(24) = {15, 18};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(15) = {17, 24, -21, -23};
Ruled Surface(10) = {15};
Line Loop(16) = {18, 23, -22, -24};
Ruled Surface(11) = {16};
			
// Defining the points
Point(19) = {-6.90988298942671, 0, 60, 10.0};
Point(20) = {0, 0, 60,10.0};
Point(21) = {6.90988298942671,0,60,10.0};
// Defining the circles
Circle(25) = {19, 20, 21};
Circle(26) = {21, 20, 19};
// Defining the surface
Line Loop(17) = {26, 25};


// Defining the lines between circle 25 and 24
Line(27) = {16, 19};
Line(28) = {18, 21};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(18) = {21, 28, -25, -27};
Ruled Surface(12) = {18};
Line Loop(19) = {22, 27, -26, -28};
Ruled Surface(13) = {19};
			
// Defining the points
Point(22) = {-7.356132180109627, 0, 70, 10.0};
Point(23) = {0, 0, 70,10.0};
Point(24) = {7.356132180109627,0,70,10.0};
// Defining the circles
Circle(29) = {22, 23, 24};
Circle(30) = {24, 23, 22};
// Defining the surface
Line Loop(20) = {30, 29};


// Defining the lines between circle 29 and 28
Line(31) = {19, 22};
Line(32) = {21, 24};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(21) = {25, 32, -29, -31};
Ruled Surface(14) = {21};
Line Loop(22) = {26, 31, -30, -32};
Ruled Surface(15) = {22};
			
// Defining the points
Point(25) = {-9.440697438826296, 0, 80, 10.0};
Point(26) = {0, 0, 80,10.0};
Point(27) = {9.440697438826296,0,80,10.0};
// Defining the circles
Circle(33) = {25, 26, 27};
Circle(34) = {27, 26, 25};
// Defining the surface
Line Loop(23) = {34, 33};


// Defining the lines between circle 33 and 32
Line(35) = {22, 25};
Line(36) = {24, 27};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(24) = {29, 36, -33, -35};
Ruled Surface(16) = {24};
Line Loop(25) = {30, 35, -34, -36};
Ruled Surface(17) = {25};
			
// Defining the points
Point(28) = {-11.968268412042981, 0, 90, 10.0};
Point(29) = {0, 0, 90,10.0};
Point(30) = {11.968268412042981,0,90,10.0};
// Defining the circles
Circle(37) = {28, 29, 30};
Circle(38) = {30, 29, 28};
// Defining the surface
Line Loop(26) = {38, 37};


// Defining the lines between circle 37 and 36
Line(39) = {25, 28};
Line(40) = {27, 30};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(27) = {33, 40, -37, -39};
Ruled Surface(18) = {27};
Line Loop(28) = {34, 39, -38, -40};
Ruled Surface(19) = {28};
			
// Defining the points
Point(31) = {-15.033297016639143, 0, 100, 10.0};
Point(32) = {0, 0, 100,10.0};
Point(33) = {15.033297016639143,0,100,10.0};
// Defining the circles
Circle(41) = {31, 32, 33};
Circle(42) = {33, 32, 31};
// Defining the surface
Line Loop(29) = {42, 41};


// Defining the lines between circle 41 and 40
Line(43) = {28, 31};
Line(44) = {30, 33};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(30) = {37, 44, -41, -43};
Ruled Surface(20) = {30};
Line Loop(31) = {38, 43, -42, -44};
Ruled Surface(21) = {31};
			
// Defining the points
Point(34) = {-17.205469890442554, 0, 110, 10.0};
Point(35) = {0, 0, 110,10.0};
Point(36) = {17.205469890442554,0,110,10.0};
// Defining the circles
Circle(45) = {34, 35, 36};
Circle(46) = {36, 35, 34};
// Defining the surface
Line Loop(32) = {46, 45};


// Defining the lines between circle 45 and 44
Line(47) = {31, 34};
Line(48) = {33, 36};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(33) = {41, 48, -45, -47};
Ruled Surface(22) = {33};
Line Loop(34) = {42, 47, -46, -48};
Ruled Surface(23) = {34};
			
// Defining the points
Point(37) = {-20.72964896828013, 0, 120, 10.0};
Point(38) = {0, 0, 120,10.0};
Point(39) = {20.72964896828013,0,120,10.0};
// Defining the circles
Circle(49) = {37, 38, 39};
Circle(50) = {39, 38, 37};
// Defining the surface
Line Loop(35) = {50, 49};


// Defining the lines between circle 49 and 48
Line(51) = {34, 37};
Line(52) = {36, 39};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(36) = {45, 52, -49, -51};
Ruled Surface(24) = {36};
Line Loop(37) = {46, 51, -50, -52};
Ruled Surface(25) = {37};
			
// Defining the points
Point(40) = {-22.212166116452387, 0, 130, 10.0};
Point(41) = {0, 0, 130,10.0};
Point(42) = {22.212166116452387,0,130,10.0};
// Defining the circles
Circle(53) = {40, 41, 42};
Circle(54) = {42, 41, 40};
// Defining the surface
Line Loop(38) = {54, 53};


// Defining the lines between circle 53 and 52
Line(55) = {37, 40};
Line(56) = {39, 42};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(39) = {49, 56, -53, -55};
Ruled Surface(26) = {39};
Line Loop(40) = {50, 55, -54, -56};
Ruled Surface(27) = {40};
			
// Defining the points
Point(43) = {-15.756957549709801, 0, 140, 10.0};
Point(44) = {0, 0, 140,10.0};
Point(45) = {15.756957549709801,0,140,10.0};
// Defining the circles
Circle(57) = {43, 44, 45};
Circle(58) = {45, 44, 43};
// Defining the surface
Line Loop(41) = {58, 57};


// Defining the lines between circle 57 and 56
Line(59) = {40, 43};
Line(60) = {42, 45};

// Defining the surfaces of the inner side of the vocal tract
Line Loop(42) = {53, 60, -57, -59};
Ruled Surface(28) = {42};
Line Loop(43) = {54, 59, -58, -60};
Ruled Surface(29) = {43};
			
// Definig the surface area of the glotis
Line Loop(44) = {58, 57};
Plane Surface(30) = {44};
