# Entrada de dados para an�lise de temperatura pelo
# m�todo dos elementos de contorno


# Matriz para defini��o de pontos que definem a geometria
# PONTOS = [n�mero do ponto, coord. x do ponto, coord. y do ponto]
function dad_2()
# Entrada de dados 


# Matriz para definiï¿½ï¿½o de pontos
# PONTO = [nï¿½mero do ponto, coordenada x do ponto, coordenada y do ponto];
L=2;
PONTOS  = [1   0   0 
          2   2*L   0 
          3   2*L  L
			 4   0  L];
          
# Segmentos que definem a geometria
#  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
#                                                  Radio] 
# Raio do segmento: > 0 -> O centro ï¿½ ï¿½ esquerda do segmento (do ponto
#                          inicial para o ponto final) 
#                   < 0 -> O centro ï¿½ ï¿½ direita do segmento (do ponto
#                          inicial para o ponto final)
#                   = 0 -> O segmento ï¿½ uma linha reta

SEGMENTOS = [1 1 2 0;
	 	   2 2 3 0;
           3 3 4 0;
	       4 4 1 0];


# Condiï¿½ï¿½es de contorno nos segmentos
# CCSeg=[Segmento,tipo da CDC em x, valor da CDC em x , ...
#                            tipo da CDC em y, valor da CDC em y]
# tipo da CDC = 0 => o deslocamento ï¿½ conhecido
# tipo da CDC = 1 => a forï¿½a de superfï¿½cie ï¿½ conhecida
# Para condiï¿½ï¿½es de contorno de forï¿½a normal conhecida proceder:
  # tipo da CDC em x = 2, valor da CDC em x = valor da forï¿½a normal
  # Neste caso pode-se atribuir quaisquer valores para tipo da CDC em y e
  # para valor da CDC em y

CCSeg=[1 1 0
    2 0 1
    3 1 0
    4 0 0];

 k=1;
 return PONTOS,SEGMENTOS,CCSeg,k

end
function dad_1(ne=10,et=1)
PONTOS  = [1   0  0
          2    1  0]

# Segmentos que definem a geometria
#  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
#                                                  Raio, tipo do elemento]
# Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
#                          inicial para o ponto final)
#                   < 0 -> O centro � � direita do segmento (do ponto
#                          inicial para o ponto final)
#                   = 0 -> O segmento � uma linha reta
# Tipo do elemento = 1 -> Elemento quadrático contínuo
#                  = 2 -> Elemento quadrático descontínuo
#                  = 3 -> Elemento linear contínuo


SEGMENTOS = [1 1 2  .5 et
         2 2 1  .5 et]

CCSeg=[1 1 0
    2 0 1];

 k=1;

	

return PONTOS,SEGMENTOS,CCSeg,k
end
function dad_0(ne=10)
PONTOS  = [1 0 0
    2 1 0
    3 1 1
    4 0 1 ]

# Segmentos que definem a geometria
#  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final
#                                                  Raio, tipo do elemento]
# Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
#                          inicial para o ponto final)
#                   < 0 -> O centro � � direita do segmento (do ponto
#                          inicial para o ponto final)
#                   = 0 -> O segmento � uma linha reta
# Tipo do elemento = 1 -> Elemento quadrático contínuo
#                  = 2 -> Elemento quadrático descontínuo


SEGMENTOS = [1 1 2 0
    2 2 3 0
    3 3 4 0
    4 4 1 0]
CCSeg=[1 1 0
    2 0 1
    3 1 0
    4 0 0];

 k=1;

return PONTOS,SEGMENTOS,CCSeg,k
end

function dad_3()

a=1;
b=2;

PONTOS  = [1 -b 0 ;
    2 b 0 ;
    3 -a 0 ;
    4 a 0];

# Segmentos que definem a geometria
#  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
#                                                  Raio]
# Raio do segmento: > 0 -> O centro é à esquerda do segmento (do ponto
#                          inicial para o ponto final)
#                   < 0 -> O centro é à direita do segmento (do ponto
#                          inicial para o ponto final)
#                   = 0 -> O segmento é uma linha reta
SEGMENTOS = [1 1 2 b;
    2 2 1 b;
    3 3 4 -a
    4 4 3 -a];

CCSeg=[1 1 -200
    2 1 -200
    3 0 100
    4 0 100];

kmat=1;

# Geração automática de pontos internos
npi=17;
NPX=npi; # Número de pontos internos na direção X
NPY=npi; # Número de pontos internos na direção Y

# disp('Erros em relação a solução analítica')
# nnos=size(NOS,1);
# npint=size(PONTOS_int,1);
# r=sqrt(PONTOS_int(:,2).^2+PONTOS_int(:,3).^2);
a=1;
b=2;
Ti=100;
qo=-200;
qi=-qo*b/a;
To=Ti-qo*b*log(b/a);

return PONTOS,SEGMENTOS,CCSeg,kmat
end
