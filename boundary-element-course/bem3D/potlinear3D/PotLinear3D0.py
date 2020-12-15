# -*- coding: utf-8 -*-
"""
Universidade de Brasília
Departamento de Engenharia Mecânica
Brasília, setembro de 2019

Programa de elementos de contorno aplicado a problemas de condução de
calor tri-dimensional sem fontes de calor concentradas

Tipo de elementos: Triangulares lineares contínuos

Autores: Éder Lima de Albuquerque
         Gustavo Silva Vaz Gontijo (ggontijo@gmail.com)

Última modificação: 30/09/2019 - 09h02min
"""

#%% BIBLIOTECAS E ARQUIVOS NECESSÁRIOS

import meshio
import entrada_de_dados
import contorno
import sistema
import numpy as np
import time
import integracao as integ


#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

t_inicio = time.time()
print('\nPrograma iniciado.')

# Lê o arquivo de entrada de dados
arquivo, CCSup, k = entrada_de_dados.dad1()
#arquivo, CCSup, k = entrada_de_dados.dad2()
#arquivo, CCSup, k = entrada_de_dados.dad3()

# Cria a malha
malha = meshio.read(arquivo + '.msh') # Lê a malha do arquivo .msh

# Cria a matriz NOS a partir da malha
NOS = malha.points
  # NOS: Matriz [NNx3] que contém as coordenadas dos vértices da malha criada,
  #      onde NN é o número de nós do problema

# Cria a matriz ELEM a partir da malha
ELEM = malha.cells_dict['triangle']
  # ELEM: Matriz [NEx3] que contém os números dos nós que formam cada elemento,
  #       onde NE é o número de elementos do problema
  
print('Número de nós:',NOS.shape[0])
print('Número de elementos:',ELEM.shape[0])


# Cria a matriz de condições de contorno dos elementos
CDC = contorno.gera_elem_cdc(malha,CCSup)
  # CDC: Matriz [NEx3] que contém a condição de contorno de cada nó dos
  #      elementos

#%% MONTAGEM E SOLUÇÃO DO SISTEMA

# Calcula as matrizes H e G

# Calcula as matrizes H e G
npg_s = 8   # Número de pontos de Gauss para a integração singular
npg_r = 6   # Número de pontos de Gauss para a integração regular
qsi,w = np.polynomial.legendre.leggauss(npg_r); # Pontos e pesos de Gauss para a integração regular (triângulo)
qsi_quad,w_quad = np.polynomial.legendre.leggauss(npg_s) # Pontos e pesos de Gauss para a integração singular (quadrilátero)
nelem = ELEM.shape[0]           # Número de elementos
nnos = NOS.shape[0]             # Número de nós
H=np.zeros((nnos, nnos))        
G=np.zeros((nnos,3*nelem))    

t_gera_dados=time.time()-t_inicio

t_matriz_inicio = time.time()
print('Calculando as matrizes H e G.')

ncores = 1      #número de threads
H[:,:],G[:,:] = integ.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad,ncores) # Cython

t_matriz=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
  # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
  # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno
print('Aplicando as condições de contorno.')
t_aplica_cdc_inicio = time.time()

A, b, T_pr = sistema.aplica_cdc(G, H, NOS, ELEM, CDC)
t_aplica_cdc = time.time()-t_aplica_cdc_inicio 

  # A: Matriz [NNxNN] contendo colunas de H e G
  # b: Vetor [NNx1] resultante da multiplicação de N colunas de H e G pelas
  #    CDC's conhecidas
print('Resolvendo o sistema linear.')

t_sistema_inicio = time.time()

# Resolve o sistema de equações A.x = b
x = np.linalg.solve(A, b)
  # x: Vetor [NNx1] que contém os termos calculados (antes desconhecidos) de
  #    temperatura e fluxo

t_sistema=time.time()-t_sistema_inicio # tempo para montagens das matrizes H e G


t_ordena_inicio=time.time()
# Separa os valores de temperatura e fluxo
print('Separando as variáveis.')

T, q = sistema.monta_Teq(NOS, ELEM, CDC, x, T_pr)
  # T: Vetor [NNx1] que contém os valores de temperatura calculados
  # q: Vetor [3NEx1] que contém os valores de fluxo calculados
t_ordena=time.time()-t_ordena_inicio

print('Gerando o arquivo de pós-processamento.')

t_saida_inicio=time.time()

# Calcula os valores de temperatura no centróide do elemento
T_centroide = contorno.calcTcentroide(T,ELEM,NOS)

# Mostra os valores de temperatura no mapa de cor
malha.cell_data_dict['gmsh:physical'][ 'Temp_med']=T_centroide
malha.point_data={"Temperatura":T}
meshio.write(arquivo + '.vtk', malha)

t_saida=time.time()-t_saida_inicio


print('Tempo para ler e gerar os dados iniciais:',t_gera_dados)
print('Tempo para montagem das matrizes H e G:',t_matriz)
print('Tempo para aplicas as condições de contorno:',t_aplica_cdc)
print('Tempo para resolver o sistema linear:',t_sistema)
print('Tempo para ordenar os dados:',t_ordena)
print('Tempo para gerar o arquivo de pós-processamento:',t_saida)
print('Tempo de processamento:',time.time()-t_inicio,'s.\nPrograma finalizado.')
