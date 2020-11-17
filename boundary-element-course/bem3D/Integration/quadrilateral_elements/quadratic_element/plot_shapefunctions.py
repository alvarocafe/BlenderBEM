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
xi=np.linspace(-1,1,npoints)
eta=xi;
Nplot=np.zeros((npoints,npoints))
XI=np.zeros((npoints,npoints))
ETA=np.zeros((npoints,npoints))
for i in range(npoints):
    for j in range(npoints):
        N=sp.compute_shapefun(xi[i],eta[j])
        XI[i,j]=xi[i]
        ETA[i,j]=eta[j]
        Nplot[i,j]=N[4]

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XI, ETA, Nplot, cmap=cm.coolwarm,linewidth=0, antialiased=False)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel(r'$\xi$')
ax.set_ylabel(r'$\eta$')
ax.set_zlabel(r'$N$')
plt.show()
