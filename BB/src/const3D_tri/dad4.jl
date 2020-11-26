function dad4()
# Programa: dad4.jl
#
# 8 elementos por face
# Todos os elementos tem a mesma orienta��o.
#
# Elementos de contorno triangulares lineares cont�nuos
#
# Cubo de lado = LL

LL = 1; # Comprimento do lado do cubo

# Coordenada dos n�s que definem a geometria (n�s geom�tricos)
# NOS = [n�mero do n�, coord. x, coord. y, coord. z];
NOS_GEO = [ 1   0  	   0  	 0
        2   LL/2   0  	 0
        3   LL 	   0	 0
        4   0      LL/2	 0
        5   LL/2   LL/2	 0
        6   LL	   LL/2	 0
        7   0	   LL    0
        8   LL/2   LL    0
        9   LL     LL    0

       10   0	   0	 LL/2
       11   LL/2   0	 LL/2
       12   LL     0	 LL/2
       13   LL     LL/2	 LL/2
       14   LL     LL    LL/2
       15   LL/2   LL    LL/2
       16   0	   LL    LL/2
       17   0	   LL/2	 LL/2
       
       18   0	   0	 LL
       19   LL/2   0	 LL
       20   LL     0	 LL
       21   0	   LL/2	 LL
       22   LL/2   LL/2	 LL
       23   LL     LL/2	 LL
       24   0	   LL    LL
       25   LL/2   LL    LL
       26   LL     LL    LL];

# Matriz de conectividade (n�s que definem os elementos)
# ELEM = [n�mero do elemento, no1, no2, no3, face] 
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

# Matriz de condi��es de contorno das faces
# CCFace = [n�mero da face, tipo da CDC, valor da CDC]
# tipo da CDC = 0 => a temperatura � conhecida
# tipo da CDC = 1 => o fluxo � conhecido
CCFace = [1 0 0
          2 1 0
          3 1 0
          4 1 0
          5 1 0
          6 0 LL];

# Fontes concentradas
# fc = [Intensidade da fonte, Coord. X, Coord. Y, Coord. Z]
# a matriz fc suporta m�ltiplas fontes concentradas; uma por linha da
# matriz
fc = [0 LL/2 LL/2 LL/2];

k = 1.; # Condutividade t�rmica do material
return NOS_GEO, ELEM, CCFace, fc, k
end
