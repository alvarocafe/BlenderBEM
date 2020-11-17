# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:38:45 2016

@author: mansour
"""
cimport numpy as np
import numpy as np 

def mount_matrix(np.ndarray[double, ndim=2] node_med,np.ndarray[double, ndim=2] normal,
                  np.ndarray[double, ndim=1] Jac,np.ndarray[double, ndim=2] nodes_coord,
                  np.ndarray[long, ndim=2] elem,double k,np.ndarray[double, ndim=2] bcs):
    cdef int nn =elem.shape[0]
    cdef np.ndarray A = np.zeros((nn, nn))
    cdef int jj,ii,no1,no2,no3,p1,p2
    cdef double X1,Y1,Z1,X2,Y2,Z2,l12,l23,m12,m23,n12,n23,XS,YS,ZS,d,a,zS,hel,gel,t12,d12,t1,t2,a1,a2,pi
    cdef np.ndarray[double, ndim=1]  vet1,vet2,vet3,o
    cdef np.ndarray b = np.zeros(nn)
    cdef np.ndarray[double, ndim=2]  nodes_e
    cdef np.ndarray[long, ndim=2]  arestas
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
        vet1=np.array([l12,m12,n12])/np.linalg.norm(np.array([l12,m12,n12]))
        vet3=np.cross(vet1,np.array([l23,m23,n23]))
        vet3=vet3/np.linalg.norm(vet3) # Vetor normal unitário ao elemento
        vet2=np.cross(vet3,vet1)
        for ii in range(0,nn):
            XS=node_med[ii,0]
            YS=node_med[ii,1]
            ZS=node_med[ii,2]
            d=np.dot(np.array([X1,Y1,Z1])-np.array([XS,YS,ZS]),vet3)/(np.linalg.norm(vet3)**2)
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
            # Dist�ncia do ponto O at� os v�rtices
                t1=l12*x1+m12*y1
                t2=l12*x2+m12*y2
            # Limites de integra��o (ap�ndice B do artigo japon�s)
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
            # Integra��o singular
                hel=-1./2.
            if bcs[jj,0] == 0:
                A[ii,jj] = -gel # Os elementos de G v�o para a matriz A
                b[ii] = b[ii] - hel*bcs[jj,1]# Os elementos de H v�o para o vetor b
            else: # O fluxo � conhecido
                A[ii,jj] = + hel # Os elementos de H v�o para a matriz A
                b[ii] = b[ii] + gel*bcs[jj,1]# Os elementos de G v�o para o vetor b
    return A,b 
