#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 08:58:56 2019

@author: eder
"""


import meshio
import entrada_de_dados


#%% ENTRADA DE DADOS E PRÉ-PROCESSAMENTO

print('\nPrograma iniciado.')

# Lê o arquivo de entrada de dados
arquivo, CCSup, k = entrada_de_dados.dad2()

# Cria a malha
malha = meshio.read(arquivo + '.msh', "gmsh4-ascii") # Lê a malha do arquivo .msh

# Cria a matriz NOS a partir da malha
NOS = malha.points
  # NOS: Matriz [NNx3] que contém as coordenadas dos vértices da malha criada,
  #      onde NN é o número de nós do problema

# Cria a matriz ELEM a partir da malha
ELEM = malha.cells['triangle']
  # ELEM: Matriz [NEx3] que contém os números dos nós que formam cada elemento,
  #       onde NE é o número de elementos do problema
  
# Cria a matriz de condições de contorno dos elementos
superf = malha.cell_data['triangle'][ 'gmsh:physical']
  
print('Número de nós:',NOS.shape[0])
print('Número de elementos:',ELEM.shape[0])

#
malha.cell_data['triangle'][ 'CDC']=superf.astype(np.float)
meshio.write(arquivo + '_cdc.vtk', malha)