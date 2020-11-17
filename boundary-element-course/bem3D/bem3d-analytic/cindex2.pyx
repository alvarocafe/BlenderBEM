# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:38:45 2016

@author: mansour
"""
cimport numpy as np
import numpy as np
from libc.math cimport log 
from libc.math cimport atan2 
from libc.math cimport sqrt 
from libc.math cimport asinh
from libc.math cimport abs

import cython

@cython.boundscheck(False)
@cython.wraparound(False)

def mount_matrix(double[:,:] node_med,double[:,:]  normal,
                  double[:] Jac,double[:,:] nodes_coord,
                 long[:,:]  elem,double k,double[:,:]  bcs):
    cdef int nn =elem.shape[0]
    cdef double[:,:] A = np.zeros((nn, nn))
    cdef int jj,ii,no1,no2,no3,p1,p2
    cdef double X1,Y1,Z1,X2,Y2,Z2,l12,l23,m12,m23,n12,n23,XS,YS,ZS,d,a,zS,hel,gel,t12,d12,t1,t2,a1,a2,pi
    cdef double[:]  vet1,vet2,vet3,o
    cdef double[:] b = np.zeros(nn)
    cdef double[:,:]   nodes_e
    cdef long[:,:]   arestas
    vet1=np.zeros(3)
    vet2=np.zeros(3)
    vet3=np.zeros(3)
    o=np.zeros(3)
    nodes_e=np.zeros((3,3))
    arestas=np.array([[0,1],[1,2],[2,0]])
    pi=np.pi
    for jj in range(0,nn):
        no1=elem[jj,0] 
        no2=elem[jj,1]                
        no3=elem[jj,2]                
        X1=nodes_coord[no1,0]
        Y1=nodes_coord[no1,1]
        Z1=nodes_coord[no1,2]
        X2=nodes_coord[no2,0]
        Y2=nodes_coord[no2,1]
        Z2=nodes_coord[no2,2]
        X3=nodes_coord[no3,0]
        Y3=nodes_coord[no3,1]
        Z3=nodes_coord[no3,2]
        l12=(X2-X1)
        l23=(X3-X2)
        m12=(Y2-Y1)
        m23=(Y3-Y2)
        n12=(Z2-Z1)
        n23=(Z3-Z2)
        vet1[0]=l12/sqrt(l12**2+m12**2+n12**2)
        vet1[1]=m12/sqrt(l12**2+m12**2+n12**2)
        vet1[2]=n12/sqrt(l12**2+m12**2+n12**2)
        vet3[0]=(vet1[1]*n23-vet1[2]*m23)
        vet3[1]=(vet1[2]*l23-vet1[0]*n23)
        vet3[2]=(vet1[0]*m23-vet1[1]*l23)
        normavet3=sqrt(vet3[0]**2+vet3[1]**2+vet3[2]**2)
        vet3[0]=vet3[0]/normavet3# Vetor normal unitário ao elemento
        vet3[1]=vet3[1]/normavet3# Vetor normal unitário ao elemento
        vet3[2]=vet3[2]/normavet3# Vetor normal unitário ao elemento
        vet2[0]=(vet3[1]*vet1[2]-vet3[2]*vet1[1])
        vet2[1]=(vet3[2]*vet1[0]-vet3[0]*vet1[2])
        vet2[2]=(vet3[0]*vet1[1]-vet3[1]*vet1[0])
        for ii in range(0,nn):
            XS=node_med[ii,0]
            YS=node_med[ii,1]
            ZS=node_med[ii,2]
            d=((X1-XS)*vet3[0]+(Y1-YS)*vet3[1]+(Z1-ZS)*vet3[2])/sqrt(vet3[0]**2+vet3[1]**2+vet3[2]**2)
            o[0]=d*vet3[0]+XS
            o[1]=d*vet3[1]+YS
            o[2]=d*vet3[2]+ZS # Coordenadas da origem do sistema local
            nodes_e[0,0]=X1
            nodes_e[0,1]=Y1
            nodes_e[0,2]=Z1
            nodes_e[1,0]=X2
            nodes_e[1,1]=Y2
            nodes_e[1,2]=Z2
            nodes_e[2,0]=X3
            nodes_e[2,1]=Y3
            nodes_e[2,2]=Z3
            zS=abs((vet3[0]*XS+vet3[1]*YS+vet3[2]*ZS)-(vet3[0]*o[0]+vet3[1]*o[1]+vet3[2]*o[2])) # Coordenada z do ponto fonte
            hel=0
            gel=0
            for i in range(3):
                p1=arestas[i,0]
                p2=arestas[i,1]
                x1=vet1[0]*nodes_e[p1,0]+vet1[1]*nodes_e[p1,1]+vet1[2]*nodes_e[p1,2]-vet1[0]*o[0]-vet1[1]*o[1]-vet1[2]*o[2]
                y1=vet2[0]*nodes_e[p1,0]+vet2[1]*nodes_e[p1,1]+vet2[2]*nodes_e[p1,2]-vet2[0]*o[0]-vet2[1]*o[1]-vet2[2]*o[2]
                x2=vet1[0]*nodes_e[p2,0]+vet1[1]*nodes_e[p2,1]+vet1[2]*nodes_e[p2,2]-vet1[0]*o[0]-vet1[1]*o[1]-vet1[2]*o[2]
                y2=vet2[0]*nodes_e[p2,0]+vet2[1]*nodes_e[p2,1]+vet2[2]*nodes_e[p2,2]-vet2[0]*o[0]-vet2[1]*o[1]-vet2[2]*o[2]
            # Comprimento da aresta
            # Cossenos diretores da aresta
                t12=np.sqrt((x2-x1)**2+(y2-y1)**2)
                l12=(x2-x1)/t12
                m12=(y2-y1)/t12
            # Dist�ncia da origem a aresta
                d12=l12*y1-m12*x1
            # Dist�ncia do ponto O at� os v�rtices
                t1=l12*x1+m12*y1
                t2=l12*x2+m12*y2
            # Limites de integra��o (ap�ndice B do artigo japon�s)
                if(d12!=0):
                    a1=t1/d12
                    a2=t2/d12        

                    hel=hel+1/(4*pi)*(-atan2(sqrt(zS**2+(a2**2+1)*d12**2),(a2*zS))+ 
                         atan2(sqrt(zS**2+(a1**2+1)*d12**2),(a1*zS))+atan2(1,a2)-atan2(1,a1))
                    gel=gel-1/(4*pi*k)*(-zS*atan2(sqrt(zS**2+(a2**2+1)*d12**2),(a2*zS))+ 
                          zS*atan2(np.sqrt(zS**2+(a1**2+1)*d12**2),(a1*zS))+ 
                          d12*asinh((a2*d12)/sqrt(zS**2+d12**2))-d12*asinh((a1*d12)/sqrt(zS**2+d12**2))+ 
                          (atan2(1,a2)-atan2(1,a1))*zS)
            if ii == jj:
            # Integra��o singular
                hel=-1./2.
            if bcs[jj,0] == 0:
                A[ii,jj] = -gel # Os elementos de G v�o para a matriz A
                b[ii] = b[ii] - hel*bcs[jj,1]# Os elementos de H v�o para o vetor b
            else: # O fluxo � conhecido
                A[ii,jj] = + hel # Os elementos de H v�o para a matriz A
                b[ii] = b[ii] + gel*bcs[jj,1]# Os elementos de G v�o para o vetor b
    return A,b 
