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
npoints=8
xiquad=np.linspace(-1,1,npoints)
a=0
b=1
xitri = 0.5*(xiquad + 1)*(b - a) + a
eta=xitri
Nplot=np.zeros((npoints,npoints))
XI=np.zeros((npoints,npoints))
ETA=np.zeros((npoints,npoints))
for i in range(npoints):
    for j in range(npoints):
        xi=(1-eta[j])*xitri[i]
        N=sp.compute_shapefun(xi,eta[j])
        XI[i,j]=xi
        ETA[i,j]=eta[j]
        Nplot[i,j]=N[0]

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XI, ETA, Nplot, cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')
ax.set_zlabel(r'$N$')
plt.show()
