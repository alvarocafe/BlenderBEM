#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 15:40:23 2017
Compute the perimeter of a cincunference using quadratic continuous elements
@author: eder
"""
import numpy as np
import matplotlib.pyplot as plt

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
    xplot=np.zeros(ngp*n_el)
    yplot=np.zeros(ngp*n_el)
    Perimeter=0
    intF=0
    for i in range(n_el):
        node1 = elem[i,0]
        node2 = elem[i,1]
        node3 = elem[i,2]
        x1 = nodes[node1,0]
        y1 = nodes[node1,1]
        x2 = nodes[node2,0]
        y2 = nodes[node2,1]
        x3 = nodes[node3,0]
        y3 = nodes[node3,1]        
        for k in range(ngp):
            N1 = 0.5*xi[k]*(xi[k]-1) # shape function N1
            N2 = 1-xi[k]**2
            N3 = 0.5*xi[k]*(xi[k]+1);
            x=N1*x1+N2*x2+N3*x3
            y=N1*y1+N2*y2+N3*y3
            ivet=ngp*i+k # index of vectors that contain points to be plotted
            xplot[ivet]=x
            yplot[ivet]=y
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

nelem = 10  # Number of elements
radius=1. # Radius of the circunferemce
f='1'

angles = np.linspace(0, 2*np.pi, 2*nelem, endpoint=False)

# Create the x and y coordinates of the points.
x = radius*np.cos(angles)
y = radius*np.sin(angles)
nodes=np.concatenate((x, y)).reshape(2,-1).transpose()
inodes=np.arange(2*nelem,dtype=int)
elem=np.zeros((nelem,3),dtype=int)
# Create the elem matrix
elem[:,0]=inodes[0::2];
elem[:,1]=inodes[1::2];
elem[:-1,2]=inodes[2::2]
elem[-1,2]=inodes[0] # last node of second column equal to first node
intf,perimeter,xplot,yplot=compute_perimeter(elem,nodes,f)
print('Perimeter = ',perimeter)
print('int f dxdy = ',intf)
xnode=np.append(nodes[:,0],nodes[0,0])
ynode=np.append(nodes[:,1],nodes[0,1])
plt.plot(xplot,yplot,"b-",markersize=4)	# Plot the node of the elements
plt.plot(xnode[::2],ynode[::2],"ro",markersize=4,label='End nodes')	# Plot the node of the elements
plt.plot(xnode[1::2],ynode[1::2],"gd",markersize=4,label='Mid nodes')	# Plot the node of the elements
plt.axis("equal")
# Create a legend of all the existing plots using their labels as names
plt.legend(loc="upper right",fancybox="true") 
