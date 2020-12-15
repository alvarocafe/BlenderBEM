# -*- coding: utf-8 -*-
"""PotLinear3D_teste.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ZtGwttcYnix9IswvaeOSaMFgnjbLYN9R
"""

from google.colab import drive
drive.mount('/content/drive')

!pip install meshio

cd /content/drive/'My Drive'/'UNB/11º Semestre'/PG2/PotLinear3D_/PotLinear3D_meshio

# In order to compile pyx files, in the terminal use the command:
!python setup.py build_ext -i
!python setup2.py build_ext -i
!python setup3.py build_ext -i
!python setup4.py build_ext -i

import contorno
import sistema
import numpy as np
import integracao as integ
import integracao2 as integ2
import integracao3 as integ3
import integracao4 as integ4

NOS =np.array([[ 0.  , 0.,   0.],
   [1. ,  0. ,  0.],
   [0.,  1 ,  0],
   [1 ,  1 ,  0],
   [0 ,  0 ,  1],
   [1  , 0 ,  1],
   [0 ,  1 ,  1],
   [1  , 1  , 1]])
# A matriz ELEM tem 4 colunas e o número de linhas é igual ao número de
#   elementos, ou seja de triângulos. Neste caso são 12 elementow
# ELEM = [ número do elemento, nó 1, nó2, nó 3]
ELEM =np.array([[1  ,  4  ,  2],
   [1 ,   3   , 4],
   [1 ,   6  ,  5],
   [1 ,   2  ,  6],
   [2 ,   8 ,   6],
   [2 ,   4 ,   8],
   [3 ,   8 ,   4],
   [3  ,  7  ,  8],
   [1  ,  7 ,   3],
   [1  ,  5  ,  7],
   [5 ,   8  ,  7],
   [5  ,  6  ,  8]])-1

CDC=np.array([[0., 0., 0., 0.],
       [0., 0., 0., 0.],
       [1., 0., 0., 0.],
       [1., 0., 0., 0.],
       [1., 0., 0., 0.],
       [1., 0., 0., 0.],
       [1., 0., 0., 0.],
       [1., 0., 0., 0.],
       [1., 0., 0., 0.],
       [1., 0., 0., 0.],
       [0., 1., 1., 1.],
       [0., 1., 1., 1.]])
k=1.
print('Número de nós:',NOS.shape[0])
print('Número de elementos:',ELEM.shape[0])
npg_s = 8   # Número de pontos de Gauss para a integração singular
npg_r = 6   # Número de pontos de Gauss para a integração regular
qsi,w = np.polynomial.legendre.leggauss(npg_r) # Pontos e pesos de Gauss para a integração regular (triângulo)
qsi_quad,w_quad = np.polynomial.legendre.leggauss(npg_s) # Pontos e pesos de Gauss para a integração singular (quadrilátero)
nelem = ELEM.shape[0]           # Número de elementos
nnos = NOS.shape[0]
H=np.zeros((nnos, nnos))        
G=np.zeros((nnos,3*nelem))    
ncores = 2                      # Número de núcleos para o integ2
# python setup.py build_ext -i

H[:,:],G[:,:] = integ.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad) # Cython
#H[:,:],G[:,:] = integ2.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad, ncores) # Cython
#H[:,:],G[:,:] = integ3.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad) # Cython
#H[:,:],G[:,:] = integ4.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad, ncores) # Cython

A, b, T_pr = sistema.aplica_cdc(G, H, NOS, ELEM, CDC)
x = np.linalg.solve(A, b)
T, q = sistema.monta_Teq(NOS, ELEM, CDC, x, T_pr)
print("x = ",x)
print("T: ",T)
print("q: ",q)
