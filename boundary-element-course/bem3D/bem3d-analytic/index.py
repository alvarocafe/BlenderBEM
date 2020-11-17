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
                hel=-1/2
            if bcs[jj,0] == 0:
                A[ii,jj] = -gel # Os elementos de G v�o para a matriz A
                b[ii] = b[ii] - hel*bcs[jj,1]# Os elementos de H v�o para o vetor b
            else: # O fluxo � conhecido
                A[ii,jj] = + hel # Os elementos de H v�o para a matriz A
                b[ii] = b[ii] + gel*bcs[jj,1]# Os elementos de G v�o para o vetor b
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

def compute_gsing(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,xiquad,wquad,k):
    npg=len(xiquad) # N�mero de pontos de integra��o
    g = 0.0 # inicializa��o da matriz G
    for kk in range(0,3):
        x1t=xd # coordenada x do primeiro n� do quadrilatero desgenerado
        y1t=yd # coordenada y do primeiro n� do quadrilatero desgenerado
        z1t=zd # coordenada z do primeiro n� do quadrilatero desgenerado
        x2t=xd # coordenada x do segundo n� do quadrilatero desgenerado
        y2t=yd # coordenada y do segundo n� do quadrilatero desgenerado
        z2t=zd # coordenada z do segundo n� do quadrilatero desgenerado
        if(kk==0): # Terceiro e quarto n�s do primeiro quadril�tero desgenerado
            x3t=x1 # coordenada x do terceiro n� do quadrilatero desgenerado
            y3t=y1 # coordenada y do terceiro n� do quadrilatero desgenerado
            z3t=z1 # coordenada z do terceiro n� do quadrilatero desgenerado
            x4t=x2 # coordenada x do quarto n� do quadrilatero desgenerado
            y4t=y2 # coordenada y do quarto n� do quadrilatero desgenerado
            z4t=z2 # coordenada z do quarto n� do quadrilatero desgenerado
        elif(kk==1): # Terceiro e quarto n�s do segundo quadril�tero desgenerado
            x3t=x2 # coordenada x do terceiro n� do quadrilatero desgenerado
            y3t=y2 # coordenada y do terceiro n� do quadrilatero desgenerado
            z3t=z2 # coordenada z do terceiro n� do quadrilatero desgenerado
            x4t=x3 # coordenada x do quarto n� do quadrilatero desgenerado
            y4t=y3 # coordenada y do quarto n� do quadrilatero desgenerado
            z4t=z3 # coordenada z do quarto n� do quadrilatero desgenerado                
        elif(kk==2):# Terceiro e quarto n�s do terceiro quadril�tero desgenerado
            x3t=x3 # coordenada x do terceiro n� do quadrilatero desgenerado
            y3t=y3 # coordenada y do terceiro n� do quadrilatero desgenerado
            z3t=z3 # coordenada z do terceiro n� do quadrilatero desgenerado
            x4t=x1 # coordenada x do quarto n� do quadrilatero desgenerado
            y4t=y1 # coordenada y do quarto n� do quadrilatero desgenerado
            z4t=z1 # coordenada z do quarto n� do quadrilatero desgenerado
        for ii in range(0,npg): # la�o sobre a primeira vari�vel de integra��o
            xi=xiquad[ii]
            for jj in range(0,npg): # la�o sobre a segunda vari�vel de integra��o
                eta=xiquad[jj]
                N1=1/4*(1-xi)*(1-eta)
                N2=1/4*(1+xi)*(1-eta)
                N3=1/4*(1+xi)*(1+eta)
                N4=1/4*(1-xi)*(1+eta)
                x = x1t*N1+x2t*N2+x3t*N3+x4t*N4 # coordenada x do
                y = y1t*N1+y2t*N2+y3t*N3+y4t*N4 # coordenada y do
                z = z1t*N1+z2t*N2+z3t*N3+z4t*N4 # coordenada z do
                dN1dxi=-1/4*(1-eta)
                dN1deta=-1/4*(1-xi)
                dN2dxi=1/4*(1-eta)
                dN2deta=-1/4*(1+xi)
                dN3dxi=1/4*(1+eta)
                dN3deta=1/4*(1+xi)
                dN4dxi=-1/4*(1+eta)
                dN4deta=1/4*(1-xi)
                    
                dxdxi = x1t*dN1dxi+x2t*dN2dxi+x3t*dN3dxi+x4t*dN4dxi
                dydxi = y1t*dN1dxi+y2t*dN2dxi+y3t*dN3dxi+y4t*dN4dxi
                dzdxi = z1t*dN1dxi+z2t*dN2dxi+z3t*dN3dxi+z4t*dN4dxi

                dxdeta = x1t*dN1deta+x2t*dN2deta+x3t*dN3deta+x4t*dN4deta
                dydeta = y1t*dN1deta+y2t*dN2deta+y3t*dN3deta+y4t*dN4deta
                dzdeta = z1t*dN1deta+z2t*dN2deta+z3t*dN3deta+z4t*dN4deta
                g1 = dydxi*dzdeta - dzdxi*dydeta;
                g2 = dzdxi*dxdeta - dxdxi*dzdeta;
                g3 = dxdxi*dydeta - dydxi*dxdeta;
                J = np.sqrt(g1**2.0 + g2**2.0 + g3**2.0);
                    
                rx=x-xd
                ry=y-yd
                rz=z-zd
                r=np.sqrt(rx**2+ry**2+rz**2)                    
                Tast = 1.0/(4.0*k*np.pi*r)
                g = g + wquad[jj]*wquad[ii]*J*Tast
    return g
