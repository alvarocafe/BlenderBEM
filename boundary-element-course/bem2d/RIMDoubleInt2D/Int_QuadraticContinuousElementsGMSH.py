#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:40:23 2017
Compute the perimeter of a cincunference using quadratic continuous elements
@author: eder
"""
import numpy as np
import matplotlib.pyplot as plt
from mesh import Mesh
from matplotlib.collections import LineCollection

def  calc_F(f,r,teta,xi,w):
    drodxi=r/2
    ngp = len(xi)
    F = 0
    for i in range(ngp):
        ro=r/2*(xi[i]+1)
        x=ro*np.cos(teta)
        y=ro*np.sin(teta)
        F=F+eval(f)*ro*drodxi*w[i]
    return F


def compute_perimeter(elem,nodes,f):
    ngp=8
    xi, weight = np.polynomial.legendre.leggauss(ngp)        
    n_el = elem.shape[0]
    xplot=np.zeros((n_el,ngp+2))
    yplot=np.zeros((n_el,ngp+2))
    Perimeter=0
    intF=0
    for i in range(n_el):
        node1 = elem[i,0]
        node3 = elem[i,1]
        node2 = elem[i,2]
        x1 = nodes[node1,0]
        y1 = nodes[node1,1]
        x2 = nodes[node2,0]
        y2 = nodes[node2,1]
        x3 = nodes[node3,0]
        y3 = nodes[node3,1]
        xplot[i,0]=x1
        yplot[i,0]=y1
        xplot[i,ngp+1]=x3
        yplot[i,ngp+1]=y3        
        for k in range(ngp):
            N1 = 0.5*xi[k]*(xi[k]-1) # shape function N1
            N2 = 1-xi[k]**2
            N3 = 0.5*xi[k]*(xi[k]+1);
            x=N1*x1+N2*x2+N3*x3
            y=N1*y1+N2*y2+N3*y3
            xplot[i,k+1]=x
            yplot[i,k+1]=y
            dN1dxi=-1/2+ xi[k]
            dN2dxi=-2*xi[k]
            dN3dxi=1/2 + xi[k]
            dxdxi=dN1dxi*x1+dN2dxi*x2+dN3dxi*x3
            dydxi=dN1dxi*y1+dN2dxi*y2+dN3dxi*y3
            dgamadxi=np.sqrt(dxdxi**2+dydxi**2)
            sx=dxdxi/dgamadxi;
            sy=dydxi/dgamadxi;
            nx=sy;
            ny=-sx;
            r=np.sqrt(x**2+y**2)
            rx=x/r
            ry=y/r
            theta=np.arctan2(ry,rx)
            nr=nx*rx+ny*ry     
            F=calc_F(f,r,theta,xi,weight)
            intF= intF+ F*nr/r * dgamadxi * weight[k]
            Perimeter = Perimeter + dgamadxi * weight[k]
    return intF,Perimeter,xplot,yplot

def show_geometry(nodes,xplot,yplot):
    xnode=nodes[:,0]
    ynode=nodes[:,1]
    # plt.plot(xplot,yplot,"bo",markersize=4)	# Plot the node of the elements
    plt.plot(xplot[:,0],yplot[:,0],"gd",markersize=8,label='End nodes')	# Plot the node of the elements
    plt.plot(xnode,ynode,"ro",markersize=2)	# Plot the node of the elements
    #plt.plot(xnode[1::2],ynode[1::2],"gd",markersize=4,label='Mid nodes')	# Plot the node of the elements
    plt.axis("equal")
    # Create a legend of all the existing plots using their labels as names
    plt.legend(loc="upper right",fancybox="true")
    ax = plt.axes()
    n_el = elem.shape[0]
    for i in range(n_el):
        x=xplot[i,:]
        y=yplot[i,:]
        line_segments = LineCollection([list(zip(x, y))],
                                        linestyles='solid')
        ax.add_collection(line_segments)
    plt.show()
    
    
    
file = 'placa_furo2'
mymesh = Mesh() #criou objeto: instanciou a classe
    #chamei o metodo do objeto
mymesh.read_msh(file + '.msh')
nodes = mymesh.Verts
elem = mymesh.Elmts[8][1]-1

f='1'
intf,perimeter,xplot,yplot=compute_perimeter(elem,nodes,f)
print('Perimeter = ',perimeter)
print('int f dxdy = ',intf)
plt.close("all")
show_geometry(nodes,xplot,yplot)


