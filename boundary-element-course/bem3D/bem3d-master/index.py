# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:38:45 2016

@author: mansour
"""
import numpy as np
#from numba import jit
#@jit
def mount_matrix(node_med,normal,Jac,nodes_coord, elem,k):
    npgauss=4
    xiquad,wquad=np.polynomial.legendre.leggauss(npgauss)
    # Translate x values from the interval [-1, 1] to [a, b]
    a=0
    b=1
    xitri = 0.5*(xiquad + 1)*(b - a) + a
    wtri=wquad*0.5*(b - a)
    nn =elem.shape[0]
    print(nn)
    H = np.zeros((nn, nn))
    G = np.zeros((nn, nn))
    for ii in range(0,nn):
        x0=node_med[ii,0]
        y0=node_med[ii,1]
        z0=node_med[ii,2]
        for jj in range(0,nn):
            no1=elem[jj,0]
            no2=elem[jj,1]                
            no3=elem[jj,2]                
            x1=nodes_coord[no1,0]
            y1=nodes_coord[no1,1]
            z1=nodes_coord[no1,2]
            x2=nodes_coord[no2,0]
            y2=nodes_coord[no2,1]
            z2=nodes_coord[no2,2]
            x3=nodes_coord[no3,0]
            y3=nodes_coord[no3,1]
            z3=nodes_coord[no3,2]
            if ii == jj:
                G[ii, jj] = compute_gsing(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,xiquad,wquad,k)
                H[ii, jj] = -0.5
            else:
                n1=normal[jj,0]
                n2=normal[jj,1]
                n3=normal[jj,2]
                intG=0
                intH=0
                for kk in range(0,npgauss):
                    for ll in range(0,npgauss):
                        eta=xitri[ll]
                        xi=(1-eta)*xitri[kk]
                        N1 = xi
                        N2 = eta
                        N3 = 1-xi-eta;
                        x=N1*x1+N2*x2+N3*x3
                        y=N1*y1+N2*y2+N3*y3
                        z=N1*z1+N2*z2+N3*z3
                        rx=x-x0;
                        ry=y-y0;
                        rz=z-z0;
                        r=np.sqrt(rx**2+ry**2+rz**2)
                        Tast = 1.0/(4.0*k*np.pi*r);
                        qast = (rx*n1 + ry*n2 + rz*n3)/(4.0*np.pi*r**3.0)
                        intG=intG+Tast*Jac[jj]*(1-eta)*wtri[kk]*wtri[ll]
                        intH=intH+qast*Jac[jj]*(1-eta)*wtri[kk]*wtri[ll]
                H[ii, jj] = intH
                G[ii, jj] = intG
    return H, G 

def mount_linear_system(H, G, bcs):
    nn=H.shape[0];
    A = np.zeros((nn, nn))
    B = np.zeros((nn, nn))
    b = np.zeros(nn)
    for column in range(0,nn):
        if bcs[column,0]== 0: # T is known
             A[:,column] = -G[:,column]
             B[:,column] = -H[:,column]
        else: #bound[line][column][0] == 1:  q is known
             A[:,column] = H[:,column]
             B[:,column] = G[:,column]
    b=B.dot(bcs[:,1])
    return A, b
  
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
