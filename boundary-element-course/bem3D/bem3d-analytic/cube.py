#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 16:01:31 2017

@author: eder
"""
import index
import numpy as np

def node_med(coord, elem):
    # element half point - list
    nelem=elem.shape[0]
    coord_med = np.zeros((nelem,3))
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
        coord_med[i,:]=([(x1 + x2 + x3)/3, (y1 + y2 + y3)/3, (z1 + z2 + z3)/3])    
        # print('element half points', node_med[line][0])    
        v1 = [x2 - x1, y2 - y1 , z2 - z1]
        v2 = [x3 - x2, y3 - y2 , z3 - z2]
        N = np.cross(v1, v2)
        Jac[i]=np.linalg.norm(N)
        normal[i,:]=(N/Jac[i])
        i=i+1
    return coord_med, normal, Jac


def calc_elem(surf_bc, surf, elem):
    elem_bc = np.zeros((len(elem), 2))
    #all elements with fluxo 0
    elem_bc[:,0] = 1
    for bc_surf,  bc_value  in surf_bc.items():
        ix = np.where(surf == bc_surf)[0]
        elem_bc[ix] = bc_value 

    return elem_bc

coord =np.array([[ 0  , 0,   0],
   [1 ,  0 ,  0],
   [0 ,  1 ,  0],
   [1 ,  1 ,  0],
   [0 ,  0 ,  1],
   [1  , 0 ,  1],
   [0 ,  1 ,  1],
   [1  , 1  , 1]])
# A matriz ELEM tem 4 colunas e o número de linhas é igual ao número de
#   elementos, ou seja de triângulos. Neste caso são 12 elementow
# ELEM = [ número do elemento, nó 1, nó2, nó 3]
elem =np.array([[1  ,  4  ,  2],
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

surf=np.array([0,0,1,1,2,2,3,3,4,4,5,5])
surf_bc={0:[0,0],5:[0,1]}
k=1.
coord_med, normal, Jac = node_med(coord, elem)


# apply boundary condition
elem_bc = calc_elem(surf_bc, surf, elem)

A,b = index.mount_matrix(coord_med,normal,Jac,coord, elem,k,elem_bc)
x = np.linalg.solve(A, b)