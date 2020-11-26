#############################################################################################################################
function dad1()
# Programa: dad1.m

# O CÁLCULO DE FLUXO DÁ CERTO

# 2 elementos por face
#
# Elementos de contorno triangulares lineares contínuos
#
# Cubo de lado = LL

LL = 1.; # Largura do cubo

# Coordenada dos nós que definem a geometria
# NOS = [número do nó, coord. x, coord. y, coord. z];
NOS_GEO = [ 1   0    0    0
        2   LL   0 	  0
        
        3   0	 LL	  0
        4   LL	 LL	  0
        
        5   0	 0	  LL
        6   LL	 0	  LL
        
        7   0	 LL	  LL
        8   LL	 LL	  LL];

# Matriz de conectividade (nós que definem os elementos)
# ELEM = [número do elemento, no1, no2, no3, face] 
ELEM = [1   1	4	2   1 
        2	1   3   4   1
        3   1   6   5   2
        4   1   2   6   2
        5   2   8   6   3
        6   2   4   8   3
        7   3   8   4   4
        8   3   7   8   4
        9   1   7   3   5
        10  1   5   7   5
        11  5   8   7   6
        12  5   6   8   6];

# Matriz de condições de contorno das faces
# CCFace = [número da face, tipo da CDC, valor da CDC]
# tipo da CDC = 0 => a temperatura é conhecida
# tipo da CDC = 1 => o fluxo é conhecido
CCFace = [1. 0. 0.
          2. 1. 0.
          3. 1. 0.
          4. 1. 0.
          5. 1. 0.
          6. 0. LL];
      

FR = 0.;
CW=1;
return NOS_GEO, ELEM,CCFace,FR,CW
end
#############################################################################################################################
function dad2()

LL = 1; # Comprimento do lado do cubo

  NOS_GEO = [ 1   0      0     0
        2   LL/2   0     0
        3   LL     0   0
        4   0      LL/2  0
        5   LL/2   LL/2  0
        6   LL     LL/2  0
        7   0    LL    0
        8   LL/2   LL    0
        9   LL     LL    0

       10   0    0   LL/2
       11   LL/2   0   LL/2
       12   LL     0   LL/2
       13   LL     LL/2  LL/2
       14   LL     LL    LL/2
       15   LL/2   LL    LL/2
       16   0    LL    LL/2
       17   0    LL/2  LL/2
       
       18   0    0   LL
       19   LL/2   0   LL
       20   LL     0   LL
       21   0    LL/2  LL
       22   LL/2   LL/2  LL
       23   LL     LL/2  LL
       24   0    LL    LL
       25   LL/2   LL    LL
       26   LL     LL    LL];

# Matriz de conectividade (nós que definem os elementos)
# ELEM = [número do elemento, no1, no2, no3, face] 
ELEM = [ 1     1     4     2     1
         2     4     5     2     1
         3     4     7     5     1
         4     7     8     5     1
         5     2     5     3     1
         6     5     6     3     1
         7     5     8     6     1
         8     8     9     6     1
         
         9     1    11    10     2
        10     1     2    11     2
        11     2    12    11     2
        12     2     3    12     2
        13    10    19    18     2
        14    10    11    19     2
        15    11    20    19     2
        16    11    12    20     2
        
        17     3    13    12     3
        18     3     6    13     3
        19     6    14    13     3
        20     6     9    14     3
        21    12    23    20     3
        22    12    13    23     3
        23    13    26    23     3
        24    13    14    26     3
        
        25     9    15    14     4
        26     9     8    15     4
        27     8    16    15     4
        28     8     7    16     4
        29    14    25    26     4
        30    14    15    25     4
        31    15    24    25     4
        32    15    16    24     4
        
        33     7    17    16     5
        34     7     4    17     5
        35     4    10    17     5
        36     4     1    10     5
        37    16    21    24     5
        38    16    17    21     5
        39    17    18    21     5
        40    17    10    18     5
        
        41    23    26    25     6
        42    23    25    22     6
        43    20    23    22     6
        44    20    22    19     6
        45    22    25    24     6
        46    22    24    21     6
        47    19    22    21     6
        48    19    21    18     6];

# Matriz de condições de contorno das faces
# CCFace = [número da face, tipo da CDC, valor da CDC]
# tipo da CDC = 0 => a temperatura é conhecida
# tipo da CDC = 1 => o fluxo é conhecido
CCFace = [1 0 0
          2 1 0
          3 1 0
          4 1 0
          5 1 0
          6 0 LL];

FR = 600.; 
CW=600.0000;
return NOS_GEO, ELEM,CCFace,FR,CW
end


function dad3()
LL = 1;

dx = LL/3;
dy = LL/3;
dz = LL/3;

  NOS_GEO = [ 1 0.0  0.0   0.0
        2 1*dx   0.0   0.0
        3 2*dx   0.0   0.0
        4 3*dx   0.0   0.0
        5 0.0  1*dy  0.0
        6 1*dx   1*dy  0.0
        7 2*dx   1*dy  0.0
        8 3*dx   1*dy  0.0
        9 0.0  2*dy  0.0
       10 1*dx   2*dy  0.0
       11 2*dx   2*dy  0.0
       12 3*dx   2*dy  0.0
       13 0.0  3*dy  0.0
       14 1*dx   3*dy  0.0
       15 2*dx   3*dy  0.0
       16 3*dx   3*dy  0.0

       17 0.0  0.0   1*dz
       18 1*dx   0.0   1*dz
       19 2*dx   0.0   1*dz
       20 3*dx   0.0   1*dz
       21 0.0  1*dy  1*dz
       22 3*dx   1*dy  1*dz
       23 0.0  2*dy  1*dz
       24 3*dx   2*dy  1*dz
       25 0.0  3*dy  1*dz
       26 1*dx   3*dy  1*dz
       27 2*dx   3*dy  1*dz
       28 3*dx   3*dy  1*dz
       
       29 0.0  0.0   2*dz
       30 1*dx   0.0   2*dz
       31 2*dx   0.0   2*dz
       32 3*dx   0.0   2*dz
       33 0.0  1*dy  2*dz
       34 3*dx   1*dy  2*dz
       35 0.0  2*dy  2*dz
       36 3*dx   2*dy  2*dz
       37 0.0  3*dy  2*dz
       38 1*dx   3*dy  2*dz
       39 2*dx   3*dy  2*dz
       40 3*dx   3*dy  2*dz

       41 0.0  0.0   3*dz
       42 1*dx   0.0   3*dz
       43 2*dx   0.0   3*dz
       44 3*dx   0.0   3*dz
       45 0.0  1*dy  3*dz
       46 1*dx   1*dy  3*dz
       47 2*dx   1*dy  3*dz
       48 3*dx   1*dy  3*dz
       49 0.0  2*dy  3*dz
       50 1*dx   2*dy  3*dz
       51 2*dx   2*dy  3*dz
       52 3*dx   2*dy  3*dz
       53 0.0  3*dy  3*dz
       54 1*dx   3*dy  3*dz
       55 2*dx   3*dy  3*dz
       56 3*dx   3*dy  3*dz];

# Matriz de conectividade (nós que definem os elementos)
# ELEM = [número do elemento, no1, no2, no3, face]
ELEM = [ 1     5     6     1     1
         2     1     6     2     1
         3     2     6     7     1
         4     2     7     3     1
         5     3     7     8     1
         6     3     8     4     1
         7     5     9    10     1
         8     5    10     6     1
         9     6    10    11     1
        10     6    11     7     1
        11     7    11    12     1
        12     7    12     8     1
        13     9    13    14     1
        14     9    14    10     1
        15    10    14    15     1
        16    10    15    11     1
        17    11    15    16     1
        18    11    16    12     1
        
        19    42    46    41     6
        20    41    46    45     6
        21    42    43    47     6
        22    42    47    46     6
        23    43    44    48     6
        24    43    48    47     6
        25    45    46    50     6
        26    45    50    49     6
        27    46    47    51     6
        28    46    51    50     6
        29    47    48    52     6
        30    47    52    51     6
        31    49    50    54     6
        32    49    54    53     6
        33    50    51    55     6
        34    50    55    54     6
        35    51    52    56     6
        36    51    56    55     6
        
        37     1    17    21     5
        38     1    21     5     5
        39     5    21    23     5
        40     5    23     9     5
        41     9    23    25     5
        42     9    25    13     5
        43    17    29    33     5
        44    17    33    21     5
        45    21    33    35     5
        46    21    35    23     5
        47    23    35    37     5
        48    23    37    25     5
        49    29    41    45     5
        50    29    45    33     5
        51    33    45    49     5
        52    33    49    35     5
        53    35    49    53     5
        54    35    53    37     5
        
        55     1     2    18     2
        56     1    18    17     2
        57     2     3    19     2
        58     2    19    18     2
        59     3     4    20     2
        60     3    20    19     2
        61    17    18    30     2
        62    17    30    29     2
        63    18    19    31     2
        64    18    31    30     2
        65    19    20    32     2
        66    19    32    31     2
        67    29    30    42     2
        68    29    42    41     2
        69    30    31    43     2
        70    30    43    42     2
        71    31    32    44     2
        72    31    44    43     2
        
        73     4     8    22     3
        74     4    22    20     3
        75     8    12    24     3
        76     8    24    22     3
        77    12    16    28     3
        78    12    28    24     3
        79    20    22    34     3
        80    20    34    32     3
        81    22    24    36     3
        82    22    36    34     3
        83    24    28    40     3
        84    24    40    36     3
        85    32    34    48     3
        86    32    48    44     3
        87    34    36    52     3
        88    34    52    48     3
        89    36    40    56     3
        90    36    56    52     3
        
        91    13    25    26     4
        92    13    26    14     4
        93    14    26    27     4
        94    14    27    15     4
        95    15    27    28     4
        96    15    28    16     4
        97    25    37    38     4
        98    25    38    26     4
        99    26    38    39     4
       100    26    39    27     4
       101    27    39    40     4
       102    27    40    28     4
       103    37    53    54     4
       104    37    54    38     4
       105    38    54    55     4
       106    38    55    39     4
       107    39    55    56     4
       108    39    56    40     4];

# Matriz de condições de contorno das faces
# CCFace = [número da face, tipo da CDC, valor da CDC]
# tipo da CDC = 0 => a temperatura é conhecida
# tipo da CDC = 1 => o fluxo é conhecido
CCFace = [1 0 0
          2 1 0
          3 1 0
          4 1 0
          5 1 0
          6 0 1];

FR = 0.;
CW=1;
return NOS_GEO, ELEM,CCFace,FR,CW
end
#############################################################################################################################
