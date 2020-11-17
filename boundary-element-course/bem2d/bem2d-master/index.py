# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:38:45 2016

@author: mansour
"""
import numpy as np

def mount_matrix(node_med,normal,nodes_coord, elem,k):
    npgauss=4
    xi, weight = np.polynomial.legendre.leggauss(npgauss)    
    nn =elem.shape[0]
    H = np.zeros((nn, nn))
    G = np.zeros((nn, nn))
    for ii in range(0,nn):
        x0=node_med[ii,0]
        y0=node_med[ii,1]
        for jj in range(0,nn):
            no1=elem[jj,0]
            no2=elem[jj,1]                
            x1=nodes_coord[no1,0]
            y1=nodes_coord[no1,1]
            x2=nodes_coord[no2,0]
            y2=nodes_coord[no2,1]
            L=np.sqrt((x2-x1)**2+(y2-y1)**2)
            if ii == jj:
                G[ii, jj] = (L/(2*np.pi*k))*(1 - np.log(L/2))
                H[ii, jj] = -0.5
            else:
                nx=normal[jj,0]
                ny=normal[jj,1]
                intG=0
                intH=0
                for kk in range(0,npgauss):
                    N1=1/2*(1-xi[kk])
                    N2=1/2*(1+xi[kk])
                    x=N1*x1+N2*x2
                    y=N1*y1+N2*y2
                    rx=x-x0;
                    ry=y-y0;
                    r=np.sqrt(rx**2+ry**2)
                    Tast=-1/(2*np.pi*k)*np.log(r)
                    qast=1/(2*np.pi*r**2)*(rx*nx+ry*ny)
                    intG=intG+Tast*L/2*weight[kk]
                    intH=intH+qast*L/2*weight[kk]
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

def mount_bcs(segmentos,btype,bc):
    nn=len(segmentos)
    bcs=np.zeros([nn,2])
    for i in range(0,nn):
        iseg=segmentos[i]
        bcs[i,0]=btype[iseg]
        bcs[i,1]=bc[iseg]
    return bcs
  


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
    