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
x=np.zeros(8)
y=np.zeros(8)
x[0]=0
y[0]=0
x[1]=2
y[1]=0
x[2]=2
y[2]=3
x[3]=0
y[3]=2
x[4]=(x[0]+x[1])/2
x[5]=(x[1]+x[2])/2
x[6]=(x[2]+x[3])/2
x[7]=(x[3]+x[0])/2
y[4]=(y[0]+y[1])/2
y[5]=(y[1]+y[2])/2
y[6]=(y[2]+y[3])/2
y[7]=(y[3]+y[0])/2
plt.figure()
plt.plot([x[0],x[4],x[1],x[5],x[2],x[6],x[3],x[7],x[0]],[y[0],y[4],y[1],y[5],y[2],y[6],y[3],y[7],y[0]],'ro-')
ax = plt.axes()
for i in range(8):
    ax.text(x[i], y[i], i+1, fontsize=15)
plt.show()

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

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XI, ETA, jacobian, cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')
ax.set_zlabel(r'$jacobian$')
