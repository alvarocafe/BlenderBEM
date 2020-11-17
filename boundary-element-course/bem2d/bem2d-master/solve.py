# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 09:47:52 2016

numerical line integration using linear quadrature quadrature

@author: mansour
"""
import numpy as np
import geometry


# linear relation (x,y) > (-1,1)
def gauss_legendre(nodes, normal, nodes_bound_coord):
    
    x1, y1 = nodes_bound_coord[nodes[0]]    
    x2, y2 = nodes_bound_coord[nodes[1]]
    
    xsi, weight = np.polynomial.legendre.leggauss(4)
    N1 = (1/2)*(1-xsi); dN1_dxsi = -1./2.
    N2 = (1/2)*(1+xsi); dN2_dxsi = 1./2.
    dx_dxsi = dN1_dxsi*x1 + dN2_dxsi*x2
    dy_dxsi = dN1_dxsi*y1 + dN2_dxsi*y2
    x=N1*x1+N2*x2;
    y=N1*y1+N2*y2;
    jacobian = np.sqrt(dx_dxsi**2 + dy_dxsi**2)
    nx=dy_dxsi/jacobian
    ny=-dx_dxsi/jacobian

    return jacobian, xsi, weight,x,y,nx,ny

def int_point(node_int,normal,nodes_coord,elem,k):
    # interior point
    npgauss=4
    xi, weight = np.polynomial.legendre.leggauss(npgauss)
    nintnodes=len(node_int)
    nnodes=len(elem)
    Hin = np.zeros((nintnodes, nnodes))
    Gin = np.zeros((nintnodes, nnodes))
    for ii in range(0,nintnodes):
        inodeint=node_int[ii]
        x0=nodes_coord[inodeint,0]
        y0=nodes_coord[inodeint,1]        
        for jj in range(0,nnodes):
            intG=0
            intH=0
            no1=elem[jj,0]
            no2=elem[jj,1]                
            x1=nodes_coord[no1,0]
            y1=nodes_coord[no1,1]
            x2=nodes_coord[no2,0]
            y2=nodes_coord[no2,1]
            L=np.sqrt((x2-x1)**2+(y2-y1)**2)
            nx=normal[jj,0]
            ny=normal[jj,1]
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
            Hin[ii, jj] = intH
            Gin[ii, jj] = intG
    return Hin, Gin

def int_point2(node_int_coord,normal,nodes_coord,elem,k):
    # interior point
    npgauss=4
    xi, weight = np.polynomial.legendre.leggauss(npgauss)
    nintnodes=len(node_int_coord)
    nnodes=len(elem)
    Hin = np.zeros((nintnodes, nnodes))
    Gin = np.zeros((nintnodes, nnodes))
    for ii in range(0,nintnodes):
        x0=node_int_coord[ii,0]
        y0=node_int_coord[ii,1]        
        for jj in range(0,nnodes):
            intG=0
            intH=0
            no1=elem[jj,0]
            no2=elem[jj,1]                
            x1=nodes_coord[no1,0]
            y1=nodes_coord[no1,1]
            x2=nodes_coord[no2,0]
            y2=nodes_coord[no2,1]
            L=np.sqrt((x2-x1)**2+(y2-y1)**2)
            nx=normal[jj,0]
            ny=normal[jj,1]
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
            Hin[ii, jj] = intH
            Gin[ii, jj] = intG
    return Hin, Gin
