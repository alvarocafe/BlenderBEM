z#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 11:39:46 2019

@author: eder
"""
cimport cython
cimport openmp
from cython.parallel cimport prange
from cython.parallel cimport parallel



import numpy as np
from libc.math cimport log 
from libc.math cimport atan2 
from libc.math cimport sqrt
from libc.math cimport asinh
from libc.stdio cimport printf

import cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # to avoid the exception checking
@cython.nonecheck(False)


def cal_HeG(double[:,:] NOS,long[:,:] ELEM,double k,double[:] qsil,double[:] w1,double[:] qsi_quad,double[:] w_quad, int ncores): 
    cdef long nelem,nnos,npg_s,npg_r,pc
    nelem = ELEM.shape[0]           # Número de elementos
    nnos = NOS.shape[0]             # Número de nós
    cdef double[:,:] H,G,N1,N2,N3,N1q,N2q,N3q,N4q
    H=np.zeros((nnos, nnos))        
    G=np.zeros((nnos,3*nelem))        
    ## = = = = = PONTOS E PESOS DE GAUSS = = = = =
    npg_s = qsi_quad.shape[0]   # Número de pontos de Gauss para a integração singular
    npg_r = qsil.shape[0]   # Número de pontos de Gauss para a integração regular
    for i in range(npg_r):
        qsil[i] = 0.5*(qsil[i] + 1)                             # Converte os pontos de Gauss para o intervalo [0, 1]
        w1[i] = w1[i]*0.5                                       # Converte os pesos de Gauss para o intervalo [0, 1]
 #   %% = = = = = FUNÇÕES DE FORMA = = = = =
    N1q = np.zeros((npg_s,npg_s))
    N2q = np.zeros((npg_s,npg_s))
    N3q = np.zeros((npg_s,npg_s))
    N4q = np.zeros((npg_s,npg_s))

    N1 = np.zeros((npg_r,npg_r))
    N2 = np.zeros((npg_r,npg_r))
    N3 = np.zeros((npg_r,npg_r))

    for l in range(0,npg_r):
        for m in range(0,npg_r):
            N1[l,m] = (1-qsil[m])*qsil[l] # qsi escrito como qsi_linha e eta
            N2[l,m] = qsil[m]
            N3[l,m] = 1-N1[l,m]-N2[l,m]
    
    for l in range(0,npg_s):
        for m in range(0,npg_s):
            N1q[l,m] = (1./4.)*(1. - qsi_quad[l])*(1. - qsi_quad[m])
            N2q[l,m] = (1./4.)*(1. + qsi_quad[l])*(1. - qsi_quad[m])
            N3q[l,m] = (1./4.)*(1. + qsi_quad[l])*(1. + qsi_quad[m])
            N4q[l,m] = (1./4.)*(1. - qsi_quad[l])*(1. + qsi_quad[m])

        
    #%% = = = = = CÁLCULO DOS ELEMENTOS DAS MATRIZES H e G = = = = =
    with nogil:
        for pc in prange(nelem, schedule='static', num_threads=ncores):
            calcula(pc,ELEM,NOS,k, N1q,N2q,N3q,N4q,N1,N2,N3,qsil,w1,qsi_quad,w_quad,npg_s,npg_r,H,G)
    for m in range(0,nnos):
        H[m,m] = 0             
        for l in range(0,nnos):
            if l != m:
                H[m,m] = H[m,m] - H[m,l]
#
    return H,G
#
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # to avoid the exception checking
@cython.nonecheck(False)
cdef void  calcula(long pc, long[:,:] ELEM, double[:,:] NOS,double k,
                   double[:,:] N1q,double[:,:] N2q,double[:,:] N3q,
                   double[:,:] N4q,double[:,:] N1,double[:,:] N2,
                   double[:,:] N3,double[:] qsil,double[:] w1,
                   double[:] qsi_quad,double[:] w_quad,long npg_s,long npg_r,
                   double[:,:] H,double[:,:] G) nogil:
    cdef long nelem,nnos,i,l,m,noglobal,ipf
    cdef long[3] nos
    nelem = ELEM.shape[0]           # Número de elementos
    nnos = NOS.shape[0]             # Número de nós
    cdef double[3] X1,X2,X3,Xd,n,g,h
    cdef double J
    nos[0] = ELEM[pc,0]      # Nós que compõem o elemento
    nos[1] = ELEM[pc,1]      # Nós que compõem o elemento
    nos[2] = ELEM[pc,2]      # Nós que compõem o elemento
    X1[0] = NOS[nos[0],0]    # Coordenadas do nó 1 do elemento
    X1[1] = NOS[nos[0],1]    # Coordenadas do nó 1 do elemento
    X1[2] = NOS[nos[0],2]    # Coordenadas do nó 1 do elemento
    X2[0] = NOS[nos[1],0]    # Coordenadas do nó 2 do elemento
    X2[1] = NOS[nos[1],1]    # Coordenadas do nó 2 do elemento
    X2[2] = NOS[nos[1],2]    # Coordenadas do nó 2 do elemento
    X3[0] = NOS[nos[2],0]    # Coordenadas do nó 3 do elemento
    X3[1] = NOS[nos[2],1]    # Coordenadas do nó 3 do elemento
    X3[2] = NOS[nos[2],2]    # Coordenadas do nó 3 do elemento
    n[0]=0
    n[1]=0
    n[2]=0        
    J=calc_vetnormal(X1,X2,X3,n)    # Vetor unitário normal ao elemento
    for pf in range(0,nnos): # Laço sobre os pontos fonte
        Xd[0] = NOS[pf,0]         # Coordenadas do ponto fonte
        Xd[1] = NOS[pf,1]         # Coordenadas do ponto fonte
        Xd[2] = NOS[pf,2]         # Coordenadas do ponto fonte
            
        if (pf==nos[0] or pf==nos[1] or pf==nos[2]): # Integração singular (o ponto fonte pertence ao elemento)
            if pf==nos[0]:   # O ponto fonte está no nó 1
                ipf = 1               # Sinaliza que o ponto fonte está sobre o primeiro nó
            elif pf==nos[1]: # O ponto fonte está no nó 2
                ipf = 2               # Sinaliza que o ponto fonte está sobre o segundo nó
            elif pf==nos[2]: # O ponto fonte está no nó 3
                ipf = 3               # Sinaliza que o ponto fonte está sobre o terceiro nó
            g[0]=0
            g[1]=0
            g[2]=0
            calcula_Gs(X1,X2,X3,Xd,ipf,N1q,N2q,N3q,N4q,qsi_quad,w_quad,k,g)
            h[0]=0
            h[1]=0
            h[2]=0
        else: # Integração regular (o ponto fonte não pertence ao elemento)
            g[0]=0
            g[1]=0
            g[2]=0
            
            h[0]=0
            h[1]=0
            h[2]=0
            calcula_HeGns(X1,X2,X3,Xd,qsil,w1,n,k,J,N1,N2,N3,g,h)
        for nolocal in range(0,3):
            noglobal = ELEM[pc,nolocal]   # Índice da matriz global H
            H[pf,noglobal] = H[pf,noglobal] + h[nolocal]
                
        G[pf,3*pc] = g[0]
        G[pf,3*pc+1] = g[1]
        G[pf,3*pc+2] = g[2]

            
#    %% = = = = = DIAGONAL DA MATRIZ H = = = = =
#     Calcula os termos da diagonal da matriz H (consideração de corpo a temperatura constante)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # to avoid the exception checking
@cython.nonecheck(False)


cdef double calc_vetnormal(double[3] X1,double[3] X2, double[3] X3,double[3] n) nogil:
    # Function que calcula o vetor unitário normal ao elemento
    cdef double[3] v1,v2
    v1[0] = X3[0] - X2[0]           # Vetor formado pela aresta 3-2 do elemento
    v1[1] = X3[1] - X2[1]           # Vetor formado pela aresta 3-2 do elemento
    v1[2] = X3[2] - X2[2]           # Vetor formado pela aresta 3-2 do elemento
    v2[0] = X1[0] - X2[0]           # Vetor formado pela aresta 1-2 do elemento
    v2[1] = X1[1] - X2[1]           # Vetor formado pela aresta 1-2 do elemento
    v2[2] = X1[2] - X2[2]           # Vetor formado pela aresta 1-2 do elemento
    n[0]=v1[1]*v2[2]-v1[2]*v2[1]
    n[1]=v1[2]*v2[0]-v1[0]*v2[2]
    n[2]=v1[0]*v2[1]-v1[1]*v2[0]
    J=sqrt(n[0]**2+n[1]**2+n[2]**2)
    n[0]=n[0]/J
    n[1]=n[1]/J
    n[2]=n[2]/J
    return J

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # to avoid the exception checking
@cython.nonecheck(False)


cdef void calcula_HeGns(double[3] X1,double[3] X2,double[3] X3,double[3] Xd,double[:] eta,double[:] w2,double[3] n,double k,double J,double[:,:] N1,double[:,:] N2,double[:,:] N3,double[3] g,double[3] h) nogil:
   """
   Integração numérica não-singular de elementos de contorno triangulares
   lineares contínuos
   """
   cdef double[3] R,Xc
   g[0] = 0;    # Inicializa o somatório de g
   g[1] = 0;    # Inicializa o somatório de g
   g[2] = 0;    # Inicializa o somatório de g
   h[0] = 0;    # Inicializa o somatório de g
   h[1] = 0;    # Inicializa o somatório de g
   h[2] = 0;    # Inicializa o somatório de g

#   n_pint1 = qsil.shape[0]    # Nro de pontos de integração na direção qsi
   cdef long n_pint
   n_pint = eta.shape[0]     # Nro de pontos de integração na direção eta
   
   cdef long l,m
   cdef double r,Tast,qast,pi
   pi = 3.141592654

   for l in range(0,n_pint):     # Laço sobre os pontos de integração
      for m in range(0,n_pint):  # Laço sobre os pontos de integração
         Xc[0] = N1[l,m]*X1[0] + N2[l,m]*X2[0] + N3[l,m]*X3[0] # coordenadas dos pontos de integração
         Xc[1] = N1[l,m]*X1[1] + N2[l,m]*X2[1] + N3[l,m]*X3[1] # coordenadas dos pontos de integração
         Xc[2] = N1[l,m]*X1[2] + N2[l,m]*X2[2] + N3[l,m]*X3[2] # coordenadas dos pontos de integração
        
         # Solução fundamental: início
         R[0] = Xc[0] - Xd[0]
         R[1] = Xc[1] - Xd[1]
         R[2] = Xc[2] - Xd[2]
         r = sqrt(R[0]**2 + R[1]**2 + R[2]**2)
         print('compilou')
         Tast = 1.0/(4.0*k*pi*r)
         qast = (R[0]*n[0] + R[1]*n[1] + R[2]*n[2])/(4.0*pi*r**3.0)
         # Solução fundamental: fim
        
         # Integral da matriz H
         h[0] += qast*N1[l,m]*(1-eta[m])*w2[l]*w2[m]*J
         h[1] += qast*N2[l,m]*(1-eta[m])*w2[l]*w2[m]*J
         h[2] += qast*N3[l,m]*(1-eta[m])*w2[l]*w2[m]*J
         # Integral da matriz G
         g[0] += Tast*N1[l,m]*(1-eta[m])*w2[l]*w2[m]*J
         g[1] += Tast*N2[l,m]*(1-eta[m])*w2[l]*w2[m]*J
         g[2] += Tast*N3[l,m]*(1-eta[m])*w2[l]*w2[m]*J
          
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # to avoid the exception checking
@cython.nonecheck(False)


cdef void calcula_Gs(double[3] X1,double[3] X2,double[3] X3,double[3] Xd, 
                     long pf,double[:,:] N1q,double[:,:] N2q,double[:,:] N3q, 
                     double[:,:] N4q,double[:]  qsi_quad, double[:] w_quad,
                     double k,double[3] g) nogil:
    """
    Function:  calcula_Gs
    Descrição: Integração singular da matriz G. O elemento triangular é
               transformado em um elemento quadrilateral degenerado. Os dois 
               primeiros nós deste quadrilátero são coincidentes (formam um 
               só vértice do triângulo) e correspondem ao ponto onde existe a
               singularidade, ou seja, ao ponto fonte. Isto faz com que haja 
               uma concentração de pontos de integração junto à 
               singularidade, além do jacobiano ser igual a zero na 
               singularidade. No caso da matriz H, os elementos da diagonal
               são calculados pela consideração de corpo a temperatura 
               constante.
    Autor:     Gustavo Gontijo, adaptado de Éder Lima de Albuquerque
       
    Última modificação: 05/05/2014 - 19h23min
    """   
    cdef long npg,l,m
    cdef double[3] R,Xc,X1t,X2t,X3t,X4t
    
    cdef double J,r,Tast,pi
    cdef double[4] dNdqsi,dNdeta
    cdef double dxdqsi,dydqsi,dzdqsi,dxdeta,dydeta,dzdeta,g1,g2,g3
    g[0] = 0           # Inicialização da matriz g
    g[1] = 0           # Inicialização da matriz g
    g[2] = 0           # Inicialização da matriz g
    npg = qsi_quad.shape[0]   # Número de pontos de integração
    pi = 3.141592654
    # Transformação do triângulo em quadrilátero
    
    X1t[0] = Xd[0]   # coordenadas do primeiro nó do quadrilatero degenerado
    X1t[1] = Xd[1]   # coordenadas do primeiro nó do quadrilatero degenerado
    X1t[2] = Xd[2]   # coordenadas do primeiro nó do quadrilatero degenerado
    
    X4t[0] = Xd[0]   # coordenadas do quarto nó do quadrilatero degenerado
    X4t[1] = Xd[1]   # coordenadas do quarto nó do quadrilatero degenerado
    X4t[2] = Xd[2]   # coordenadas do quarto nó do quadrilatero degenerado
    
    if (pf == 1):   # O ponto fonte está no nó 1
        X2t[0] = X2[0]             # coordenadas do segundo nó do quadrilatero degenerado
        X2t[1] = X2[1]             # coordenadas do segundo nó do quadrilatero degenerado
        X2t[2] = X2[2]             # coordenadas do segundo nó do quadrilatero degenerado

        X3t[0] = X3[0]             # coordenadas do terceiro nó do quadrilatero degenerado
        X3t[1] = X3[1]             # coordenadas do terceiro nó do quadrilatero degenerado
        X3t[2] = X3[2]             # coordenadas do terceiro nó do quadrilatero degenerado
    elif (pf == 2): # O ponto fonte está no nó 2
        X2t[0] = X3[0]             # coordenadas do segundo nó do quadrilatero degenerado
        X2t[1] = X3[1]             # coordenadas do segundo nó do quadrilatero degenerado
        X2t[2] = X3[2]             # coordenadas do segundo nó do quadrilatero degenerado

        X3t[0] = X1[0]             # coordenadas do terceiro nó do quadrilatero degenerado
        X3t[1] = X1[1]             # coordenadas do terceiro nó do quadrilatero degenerado
        X3t[2] = X1[2]             # coordenadas do terceiro nó do quadrilatero degenerado
    elif (pf == 3): # O ponto fonte está no nó 3
        X2t[0] = X1[0]             # coordenada do segundo nó do quadrilatero degenerado
        X2t[1] = X1[1]             # coordenada do segundo nó do quadrilatero degenerado
        X2t[2] = X1[2]             # coordenada do segundo nó do quadrilatero degenerado

        X3t[0] = X2[0]             # coordenada do terceiro nó do quadrilatero degenerado
        X3t[1] = X2[1]             # coordenada do terceiro nó do quadrilatero degenerado
        X3t[2] = X2[2]             # coordenada do terceiro nó do quadrilatero degenerado
    
#     Integração regular do elemento
    for l in range(0,npg):       # Laço sobre a primeira variável de integração
        for m in range(0,npg):   # Laço sobre a segunda variável de integração
            # Cálculo das funções de forma
            # Coordenadas do ponto campo
            Xc[0] = N1q[l,m]*X1t[0] + N2q[l,m]*X2t[0] + N3q[l,m]*X3t[0] + N4q[l,m]*X4t[0] # coordenadas dos pontos de integração
            Xc[1] = N1q[l,m]*X1t[1] + N2q[l,m]*X2t[1] + N3q[l,m]*X3t[1] + N4q[l,m]*X4t[1] # coordenadas dos pontos de integração
            Xc[2] = N1q[l,m]*X1t[2] + N2q[l,m]*X2t[2] + N3q[l,m]*X3t[2] + N4q[l,m]*X4t[2]# coordenadas dos pontos de integração

            # Jacobiano (varia ao longo do elemento degenerado)
            dNdqsi[0] = (1./4.)*(-(1. - qsi_quad[m]))
            dNdqsi[1] = (1./4.)*(1.-qsi_quad[m])
            dNdqsi[2] = (1./4.)*(1. + qsi_quad[m])
            dNdqsi[3] = (1./4.)*(-(1. + qsi_quad[m]))
           
            dNdeta[0] = (1./4.)*(-(1. - qsi_quad[l]))
            dNdeta[1] = (1./4.)*(-(1. + qsi_quad[l]))
            dNdeta[2] = (1./4.)*((1. + qsi_quad[l]))
            dNdeta[3] = (1./4.)*((1. - qsi_quad[l]))
   
   
            dxdqsi = X1t[0]*dNdqsi[0] + X2t[0]*dNdqsi[1] + X3t[0]*dNdqsi[2] + X4t[0]*dNdqsi[3]
            dydqsi = X1t[1]*dNdqsi[0] + X2t[1]*dNdqsi[1] + X3t[1]*dNdqsi[2] + X4t[1]*dNdqsi[3]
            dzdqsi = X1t[2]*dNdqsi[0] + X2t[2]*dNdqsi[1] + X3t[2]*dNdqsi[2] + X4t[2]*dNdqsi[3]
   
            dxdeta = X1t[0]*dNdeta[0] + X2t[0]*dNdeta[1] + X3t[0]*dNdeta[2] + X4t[0]*dNdeta[3]
            dydeta = X1t[1]*dNdeta[0] + X2t[1]*dNdeta[1] + X3t[1]*dNdeta[2] + X4t[1]*dNdeta[3]
            dzdeta = X1t[2]*dNdeta[0] + X2t[2]*dNdeta[1] + X3t[2]*dNdeta[2] + X4t[2]*dNdeta[3]
           
            g1 = dydqsi*dzdeta - dzdqsi*dydeta
            g2 = dzdqsi*dxdeta - dxdqsi*dzdeta
            g3 = dxdqsi*dydeta - dydqsi*dxdeta
            J = sqrt(g1**2.0 + g2**2.0 + g3**2.0)

            # Solução fundamental: início
            R[0] = Xc[0] - Xd[0]
            R[1] = Xc[1] - Xd[1]
            R[2] = Xc[2] - Xd[2]

            r = sqrt(R[0]**2 + R[1]**2 + R[2]**2)
            Tast = 1.0/(4.0*k*pi*r)
            # Solução fundamental: fim

            # Integral da matriz G
            if pf==1:
                g[0] += Tast*(N1q[l,m]+N4q[l,m])*w_quad[l]*w_quad[m]*J
                g[1] += Tast*N2q[l,m]*w_quad[l]*w_quad[m]*J
                g[2] += Tast*N3q[l,m]*w_quad[l]*w_quad[m]*J
            elif pf==2:
                g[0] += Tast*N3q[l,m]*w_quad[l]*w_quad[m]*J
                g[1] += Tast*(N1q[l,m]+N4q[l,m])*w_quad[l]*w_quad[m]*J
                g[2] += Tast*N2q[l,m]*w_quad[l]*w_quad[m]*J
            else:
                g[0] += Tast*N2q[l,m]*w_quad[l]*w_quad[m]*J
                g[1] += Tast*N3q[l,m]*w_quad[l]*w_quad[m]*J
                g[2] += Tast*(N1q[l,m]+N4q[l,m])*w_quad[l]*w_quad[m]*J



