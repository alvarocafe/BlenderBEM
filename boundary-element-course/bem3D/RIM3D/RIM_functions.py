# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 22:38:45 2016

@author: mansour
"""
import numpy as np
from numpy import linalg as LA

def compute_normal(coord, elem):
    # element half point - list
    nelem=elem.shape[0]
    normal = np.zeros((nelem,3))
    Jac = np.zeros(nelem)
    i=0
    # loop through elements
    for node1, node2, node3 in elem:
        # first node
        x1, y1, z1 = coord[node1]
    #    print(x1,y1)
        # second node
        x2, y2, z2 = coord[node2]
    #    print(x2,y2)
        # third node
        x3, y3, z3 = coord[node3]
        # font nodes
        v1 = [x2 - x1, y2 - y1 , z2 - z1]
        v2 = [x3 - x2, y3 - y2 , z3 - z2]
        N = np.cross(v1, v2)
        Jac[i]=LA.norm(N)
        normal[i,:]=(N/Jac[i])
        i=i+1
    return normal, Jac


def compute_geoprop(normal,Jac,nodes_coord, elem,k):
    npgauss=4
    xiquad,wquad=np.polynomial.legendre.leggauss(npgauss)
    # Translate x values from the interval [-1, 1] to [a, b]
    a=0
    b=1
    xitri = 0.5*(xiquad + 1)*(b - a) + a
    wtri=wquad*0.5*(b - a)
    nn =elem.shape[0]
    surface_area = 0
    volume = 0
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
        elem_area=Jac[jj]/2
        surface_area=surface_area+elem_area
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
                rx=x
                ry=y
                rz=z
                r=np.array([rx,ry,rz])
                nr=np.dot(normal[jj,:],r)
                volume=volume+nr/3*Jac[jj]*(1-eta)*wtri[kk]*wtri[ll]
    return surface_area,volume 
