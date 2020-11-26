function dad_0()
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


SEGMENTOS = [1 1 2 0
    2 2 3 0
    3 3 4 0
    4 4 1 0]

# Matriz para defini��o da malha

# MALHA =[numero do segmento, numero de elementos no segmento]
ne=30
MALHA = [1 ne
    2 ne
    3 ne
    4 ne]

# Condi��es de contorno nos segmentos
# CCSeg=[no do segmento, tipo da CDC, valor da CDC]
# tipo da CDC = 0 => a temperatura � conhecida
# tipo da CDC = 1 => O fluxo � conhecido
CCSeg=[1 1 0
    2 0 1
    3 1 0
    4 0 0]


# Condutividade T�rmica do material
k=1
return PONTOS, SEGMENTOS, MALHA, CCSeg, k
end
