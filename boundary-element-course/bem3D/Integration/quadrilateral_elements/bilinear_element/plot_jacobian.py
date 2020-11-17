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
x=np.zeros(4)
y=np.zeros(4)
x[0]=0
y[0]=0
x[1]=2
y[1]=0
x[2]=2
y[2]=3
x[3]=0
y[3]=2
npoints=8
xi=np.linspace(-1,1,npoints)
eta=xi;
jacobian=np.zeros((npoints,npoints))
XI=np.zeros((npoints,npoints))
ETA=np.zeros((npoints,npoints))
for i in range(npoints):
    for j in range(npoints):
        dN=sp.compute_dshapefun(xi[i],eta[j])
        dxdxi=np.dot(dN[0,:],x)
        dxdeta=np.dot(dN[1,:],x)
        dydxi=np.dot(dN[0,:],y)
        dydeta=np.dot(dN[1,:],y)
        J=dxdxi*dydeta-dxdeta*dydxi;        
        XI[i,j]=xi[i]
        ETA[i,j]=eta[j]
        jacobian[i,j]=J
plt.figure()
plt.plot(np.append(x,x[0]),np.append(y,y[0]),'ro-')

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XI, ETA, jacobian, cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')
ax.set_zlabel(r'$jacobian$')
plt.show()
