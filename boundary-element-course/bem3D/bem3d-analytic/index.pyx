# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:38:45 2016

@author: mansour
"""
import numpy as np
from numpy import linalg as LA
from numba import jit

@jit
def mount_matrix(node_med,normal,Jac,nodes_coord, elem,k,bcs):
    nn =elem.shape[0]
    A = np.zeros((nn, nn))
    b = np.zeros(nn)
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
        vet1=np.array([l12,m12,n12])/LA.norm(np.array([l12,m12,n12]))
        vet3=np.cross(vet1,np.array([l23,m23,n23]))
        vet3=vet3/LA.norm(vet3) # Vetor normal unitário ao elemento
        vet2=np.cross(vet3,vet1)
        for ii in range(0,nn):
            XS=node_med[ii,0]
            YS=node_med[ii,1]
            ZS=node_med[ii,2]
            d=np.dot(np.array([X1,Y1,Z1])-np.array([XS,YS,ZS]),vet3)/(LA.norm(vet3)**2)
            o=np.dot(d,vet3)*np.ones(3)+np.array([XS,YS,ZS]) # Coordenadas da origem do sistema local
            arestas=np.array([[0,1],[1,2],[2,0]]) # N�mero dos n�s das 3 arestas do elemento
            nodes_e=np.array([[X1,Y1,Z1],[X2,Y2,Z2],[X3,Y3,Z3]])
            zS=np.abs(np.dot(vet3,[XS,YS,ZS])-np.dot(vet3,o)) # Coordenada z do ponto fonte
            hel=0
            gel=0
            for i in range(3):
                p1=arestas[i,0]
                p2=arestas[i,1]
                x1=np.dot(vet1,[nodes_e[p1,0],nodes_e[p1,1],nodes_e[p1,2]])-np.dot(vet1,o)
                y1=np.dot(vet2,[nodes_e[p1,0],nodes_e[p1,1],nodes_e[p1,2]])-np.dot(vet2,o)
                x2=np.dot(vet1,[nodes_e[p2,0],nodes_e[p2,1],nodes_e[p2,2]])-np.dot(vet1,o)
                y2=np.dot(vet2,[nodes_e[p2,0],nodes_e[p2,1],nodes_e[p2,2]])-np.dot(vet2,o)
            # Comprimento da aresta
            # Cossenos diretores da aresta
                t12=np.sqrt((x2-x1)**2+(y2-y1)**2)
                l12=(x2-x1)/t12
                m12=(y2-y1)/t12
            # Dist�ncia da origem a aresta
                d12=l12*y1-m12*x1
            # Dist�ncia do ponto O até os vértices
                t1=l12*x1+m12*y1
                t2=l12*x2+m12*y2
            # Limites de integração (apêndice B do artigo japonês)
                if(d12!=0):
                    a1=t1/d12
                    a2=t2/d12        
                    hel=hel+1/(4*np.pi)*(-np.arctan2(np.sqrt(zS**2+(a2**2+1)*d12**2),(a2*zS))+ 
                         np.arctan2(np.sqrt(zS**2+(a1**2+1)*d12**2),(a1*zS))+np.arctan2(1,a2)-np.arctan2(1,a1))
                    gel=gel-1/(4*np.pi*k)*(-zS*np.arctan2(np.sqrt(zS**2+(a2**2+1)*d12**2),(a2*zS))+ 
                          zS*np.arctan2(np.sqrt(zS**2+(a1**2+1)*d12**2),(a1*zS))+ 
                          d12*np.arcsinh((a2*d12)/np.sqrt(zS**2+d12**2))-d12*np.arcsinh((a1*d12)/np.sqrt(zS**2+d12**2))+ 
                          (np.arctan2(1,a2)-np.arctan2(1,a1))*zS)
            if ii == jj:
            # Integração singular
                hel=-1/2
            if bcs[jj,0] == 0:
                A[ii,jj] = -gel # Os elementos de G vão para a matriz A
                b[ii] = b[ii] - hel*bcs[jj,1]# Os elementos de H vão para o vetor b
            else: # O fluxo � conhecido
                A[ii,jj] = + hel # Os elementos de H vão para a matriz A
                b[ii] = b[ii] + gel*bcs[jj,1]# Os elementos de G vão para o vetor b
    return A,b 
  
def mount_vector(z, bcs):
    # mount the boundary T and q
    nn=len(z)
    T = np.zeros(nn)
    q = np.zeros(nn)
    for elem in range(0,nn):
        if bcs[elem,0] == 0:
            T[elem] = bcs[elem,1]
            q[elem] = z[elem]
        else:
            T[elem] = z[elem]
            q[elem] = bcs[elem,1]
    return T, q