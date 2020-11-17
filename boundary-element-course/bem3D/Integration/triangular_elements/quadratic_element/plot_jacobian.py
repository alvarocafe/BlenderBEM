#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 22:00:21 2017

@author: eder
"""
import numpy as np
import shape_functions as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
plt.close("all")
x=np.zeros(9)
y=np.zeros(9)
x[0]=0
y[0]=0
x[1]=2
y[1]=0
x[2]=2
y[2]=3
x[3]=0
y[3]=2.
x[4]=(x[0]+x[1])/2
x[5]=(x[1]+x[2])/2
x[6]=(x[2]+x[3])/2
x[7]=(x[3]+x[0])/2
y[4]=(y[0]+y[1])/2
y[5]=(y[1]+y[2])/2
y[6]=(y[2]+y[3])/2
y[7]=(y[3]+y[0])/2
x[8]=1
y[8]=1.5
nnos=x.shape[0]

elem=np.array([[0,1,2,4,5,8],[0,2,3,8,6,7]],dtype=int)
nelem=elem.shape[0]

plt.figure()
ax = plt.axes()
for i in range(nnos):
    ax.text(x[i], y[i], i+1, fontsize=15)
plt.show()


npoints=8
xiquad=np.linspace(-1,1,npoints)
# Transformation of Gauss quadrature to [a,b]
a=0
b=1
xitri = 0.5*(xiquad + 1)*(b - a) + a
eta=xitri
XI=np.zeros((npoints,npoints))
ETA=np.zeros((npoints,npoints))
jacobian=np.zeros((nelem,npoints,npoints))
for k in range(nelem): # loop over elements
    nos=elem[k,:]
    xnodes=x[nos]
    ynodes=y[nos]
    for i in range(npoints): # loop over eta
        for j in range(npoints): # loop over xitri
            xi=(1-eta[j])*xitri[i]
            N=sp.compute_shapefun(xi,eta[j])
            dN=sp.compute_dshapefun(xi,eta[j])
            xx=np.dot(N,xnodes)
            yy=np.dot(N,ynodes)
            dxdxi=np.dot(dN[0,:],xnodes)
            dxdeta=np.dot(dN[1,:],xnodes)
            dydxi=np.dot(dN[0,:],ynodes)
            dydeta=np.dot(dN[1,:],ynodes)
            J=dxdxi*dydeta-dxdeta*dydxi
            XI[i,j]=xi
            ETA[i,j]=eta[j]
            jacobian[k,i,j]=J
            plt.plot(np.append(xnodes[0:3],xnodes[0]),np.append(ynodes[0:3],ynodes[0]),'ro-')
            plt.plot(xnodes,ynodes,'gd')          
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XI, ETA, jacobian[1,:,:], cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')
ax.set_zlabel(r'$jacobian$')
plt.show()
