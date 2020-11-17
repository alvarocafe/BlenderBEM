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
plt.close("All")
x=np.zeros(8)
y=np.zeros(8)
z=np.zeros(8)
x[0]=0
y[0]=0
z[0]=0
x[1]=2
y[1]=0
z[1]=0
x[2]=2
y[2]=3
z[2]=0
x[3]=0
y[3]=2
z[3]=0
x[4]=(x[0]+x[1])/2
x[5]=(x[1]+x[2])/2
x[6]=(x[2]+x[3])/2
x[7]=(x[3]+x[0])/2
y[4]=(y[0]+y[1])/2
y[5]=(y[1]+y[2])/2
y[6]=(y[2]+y[3])/2
y[7]=(y[3]+y[0])/2
z[4]=(z[0]+z[1])/2
z[5]=(z[1]+z[2])/2
z[6]=(z[2]+z[3])/2
z[7]=(z[3]+z[0])/2
npoints=8
xi=np.linspace(-1,1,npoints)
eta=xi;
jacobian=np.zeros((npoints,npoints))
XI=np.zeros((npoints,npoints))
ETA=np.zeros((npoints,npoints))
X=np.zeros((npoints,npoints))
Y=np.zeros((npoints,npoints))
Z=np.zeros((npoints,npoints))
for i in range(npoints):
    for j in range(npoints):
        N=sp.compute_shapefun(xi[i],eta[j])        
        dN=sp.compute_dshapefun(xi[i],eta[j])
        xx=np.dot(N,x)
        yy=np.dot(N,y)
        zz=np.dot(N,z)
        dxdxi=np.dot(dN[0,:],x)
        dxdeta=np.dot(dN[1,:],x)
        dydxi=np.dot(dN[0,:],y)
        dydeta=np.dot(dN[1,:],y)
        dzdxi=np.dot(dN[0,:],z)
        dzdeta=np.dot(dN[1,:],z)
        
        g1=dydxi*dzdeta-dydeta*dzdxi;
        g2=dxdxi*dzdeta-dxdeta*dzdxi;
        g3=dxdxi*dydeta-dxdeta*dydxi;
        J=np.sqrt(g1**2+g2**2+g3**2)
        XI[i,j]=xi[i]
        ETA[i,j]=eta[j]
        jacobian[i,j]=J
        X[i,j]=xx
        Y[i,j]=yy
        Z[i,j]=zz

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$z$')

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XI, ETA, jacobian, cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')
ax.set_zlabel(r'$jacobian$')
