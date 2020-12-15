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
import integracao2 as integ2
import integracao3 as integ3


it=10        #número de iterações
casos=4     #número de casos com qtd. de núcleo diferente

t_matriz_1 = np.zeros((casos+1, it))
t_matriz_2 = np.zeros((casos+1, it))
t_matriz_3 = np.zeros((casos+1, it))

print (t_matriz_1)
print ('Número de iterações:', it)
	
#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

t_inicio = time.time()
print('\nPrograma iniciado.')

# Lê o arquivo de entrada de dados
arquivo, CCSup, k = entrada_de_dados.dad1()
print('Utilizando a malha 1')

#arquivo, CCSup, k = entrada_de_dados.dad2()
#print('Utilizando a malha 2')

#arquivo, CCSup, k = entrada_de_dados.dad3()
#print('Utilizando a malha 3')

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

print('Calculando código serial')
for j in range (it):
	t_gera_dados=time.time()-t_inicio

	t_matriz_inicio = time.time()

	H[:,:],G[:,:] = integ3.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad) # Cython

	t_matriz_1[0, j-1]=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
	  # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
	  # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno


for i in range (casos):
	
	ncores = i*2+2      #número de threads
	#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

	t_inicio = time.time()
	print('\nPrograma iniciado.')

	# Lê o arquivo de entrada de dados
	arquivo, CCSup, k = entrada_de_dados.dad1()
	print('Utilizando a malha 1')

	#arquivo, CCSup, k = entrada_de_dados.dad2()
	#print('Utilizando a malha 2')

	#arquivo, CCSup, k = entrada_de_dados.dad3()
	#print('Utilizando a malha 3')

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

	print('Calculando para n=',ncores)
	for j in range (it):
		t_gera_dados=time.time()-t_inicio

		t_matriz_inicio = time.time()

		H[:,:],G[:,:] = integ.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad,ncores) # Cython

		t_matriz_1[i+1, j-1]=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
		  # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
		  # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno

t_serial_med = np.mean(t_matriz_1[0,:])
t_serial_std = np.std(t_matriz_1[0,:])

t_n2_med = np.mean(t_matriz_1[1,:])
t_n2_std = np.std(t_matriz_1[1,:])

t_n4_med = np.mean(t_matriz_1[2,:])
t_n4_std = np.std(t_matriz_1[2,:])

t_n6_med = np.mean(t_matriz_1[3,:])
t_n6_std = np.std(t_matriz_1[3,:])

t_n8_med = np.mean(t_matriz_1[4,:])
t_n8_std = np.std(t_matriz_1[4,:])

print('\n tempo serial médio =', np.round(t_serial_med,2),'  desvio padrão =', np.round(t_serial_std,2))
print('\n tempo médio 2 núcleos =', np.round(t_n2_med,2),'  desvio padrão =', np.round(t_n2_std,2))
print('\n tempo médio 4 núcleos =', np.round(t_n4_med,2),'  desvio padrão =', np.round(t_n4_std,2))
print('\n tempo médio 6 núcleos =', np.round(t_n6_med,2),'  desvio padrão =', np.round(t_n6_std,2))
print('\n tempo médio 8 núcleos =', np.round(t_n8_med,2),'  desvio padrão =', np.round(t_n8_std,2))


#############################################################################################################

#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

t_inicio = time.time()
print('\nPrograma iniciado.')

# Lê o arquivo de entrada de dados
#arquivo, CCSup, k = entrada_de_dados.dad1()
#print('Utilizando a malha 1')

arquivo, CCSup, k = entrada_de_dados.dad2()
print('Utilizando a malha 2')

#arquivo, CCSup, k = entrada_de_dados.dad3()
#print('Utilizando a malha 3')

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

print('Calculando código serial')
for j in range (it):
	t_gera_dados=time.time()-t_inicio

	t_matriz_inicio = time.time()

	H[:,:],G[:,:] = integ3.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad) # Cython

	t_matriz_1[0, j-1]=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
	  # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
	  # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno


for i in range (casos):
	
	ncores = i*2+2      #número de threads
	#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

	t_inicio = time.time()
	print('\nPrograma iniciado.')

	# Lê o arquivo de entrada de dados
	#arquivo, CCSup, k = entrada_de_dados.dad1()
	#print('Utilizando a malha 1')

	arquivo, CCSup, k = entrada_de_dados.dad2()
	print('Utilizando a malha 2')

	#arquivo, CCSup, k = entrada_de_dados.dad3()
	#print('Utilizando a malha 3')

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

	print('Calculando para n=',ncores)
	for j in range (it):
		t_gera_dados=time.time()-t_inicio

		t_matriz_inicio = time.time()

		H[:,:],G[:,:] = integ.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad,ncores) # Cython

		t_matriz_1[i+1, j-1]=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
		  # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
		  # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno

t_serial_med = np.mean(t_matriz_1[0,:])
t_serial_std = np.std(t_matriz_1[0,:])

t_n2_med = np.mean(t_matriz_1[1,:])
t_n2_std = np.std(t_matriz_1[1,:])

t_n4_med = np.mean(t_matriz_1[2,:])
t_n4_std = np.std(t_matriz_1[2,:])

t_n6_med = np.mean(t_matriz_1[3,:])
t_n6_std = np.std(t_matriz_1[3,:])

t_n8_med = np.mean(t_matriz_1[4,:])
t_n8_std = np.std(t_matriz_1[4,:])

print('\n tempo serial médio =', np.round(t_serial_med,2),'  desvio padrão =', np.round(t_serial_std,2))
print('\n tempo médio 2 núcleos =', np.round(t_n2_med,2),'  desvio padrão =', np.round(t_n2_std,2))
print('\n tempo médio 4 núcleos =', np.round(t_n4_med,2),'  desvio padrão =', np.round(t_n4_std,2))
print('\n tempo médio 6 núcleos =', np.round(t_n6_med,2),'  desvio padrão =', np.round(t_n6_std,2))
print('\n tempo médio 8 núcleos =', np.round(t_n8_med,2),'  desvio padrão =', np.round(t_n8_std,2))


####################################################################################################

#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

t_inicio = time.time()
print('\nPrograma iniciado.')

# Lê o arquivo de entrada de dados
#arquivo, CCSup, k = entrada_de_dados.dad1()
#print('Utilizando a malha 1')

#arquivo, CCSup, k = entrada_de_dados.dad2()
#print('Utilizando a malha 2')

arquivo, CCSup, k = entrada_de_dados.dad3()
print('Utilizando a malha 3')

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

print('Calculando código serial')
for j in range (it):
	t_gera_dados=time.time()-t_inicio

	t_matriz_inicio = time.time()

	H[:,:],G[:,:] = integ3.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad) # Cython

	t_matriz_1[0, j-1]=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
	  # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
	  # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno


for i in range (casos):
	
	ncores = i*2+2      #número de threads
	#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

	t_inicio = time.time()
	print('\nPrograma iniciado.')

	# Lê o arquivo de entrada de dados
	#arquivo, CCSup, k = entrada_de_dados.dad1()
	#print('Utilizando a malha 1')

	#arquivo, CCSup, k = entrada_de_dados.dad2()
	#print('Utilizando a malha 2')

	arquivo, CCSup, k = entrada_de_dados.dad3()
	print('Utilizando a malha 3')

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

	print('Calculando para n=',ncores)
	for j in range (it):
		t_gera_dados=time.time()-t_inicio

		t_matriz_inicio = time.time()

		H[:,:],G[:,:] = integ.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad,ncores) # Cython

		t_matriz_1[i+1, j-1]=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
		  # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
		  # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno

t_serial_med = np.mean(t_matriz_1[0,:])
t_serial_std = np.std(t_matriz_1[0,:])

t_n2_med = np.mean(t_matriz_1[1,:])
t_n2_std = np.std(t_matriz_1[1,:])

t_n4_med = np.mean(t_matriz_1[2,:])
t_n4_std = np.std(t_matriz_1[2,:])

t_n6_med = np.mean(t_matriz_1[3,:])
t_n6_std = np.std(t_matriz_1[3,:])

t_n8_med = np.mean(t_matriz_1[4,:])
t_n8_std = np.std(t_matriz_1[4,:])

print('\n tempo serial médio =', np.round(t_serial_med,2),'  desvio padrão =', np.round(t_serial_std,2))
print('\n tempo médio 2 núcleos =', np.round(t_n2_med,2),'  desvio padrão =', np.round(t_n2_std,2))
print('\n tempo médio 4 núcleos =', np.round(t_n4_med,2),'  desvio padrão =', np.round(t_n4_std,2))
print('\n tempo médio 6 núcleos =', np.round(t_n6_med,2),'  desvio padrão =', np.round(t_n6_std,2))
print('\n tempo médio 8 núcleos =', np.round(t_n8_med,2),'  desvio padrão =', np.round(t_n8_std,2))


