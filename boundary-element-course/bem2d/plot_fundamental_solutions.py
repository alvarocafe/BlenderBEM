#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 08:28:41 2017

@author: eder
"""
'''
=================================
3D surface with polar coordinates
=================================

Demonstrates plotting a surface defined in polar coordinates.
Uses the reversed version of the YlGnBu color map.
Also demonstrates writing axis labels with latex math mode.

Example contributed by Armin Moser.
'''

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import numpy as np

plt.close("all")
# Create the mesh in polar coordinates and compute corresponding Z.
r = np.linspace(0.15, 1.25, 50)
theta = np.linspace(0, 2*np.pi, 50)
R, THETA = np.meshgrid(r, theta)
Rx=R*np.cos(THETA)
Ry=R*np.sin(THETA)
nx=np.sqrt(2)/2
ny=np.sqrt(2)/2
Tast = -1/(2*np.pi)*np.log(R)
qast=1/(2*np.pi*R**2)*(Rx*nx+Ry*ny)
# Express the mesh in the cartesian system.
X, Y = R*np.cos(THETA), R*np.sin(THETA)

# Plot the surface.

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Tast, cmap=plt.cm.hsv)
# Tweak the limits and add latex math labels.
ax.set_zlim(0, 1)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$T^*$')
plt.title("Temperature fundamental solution")
plt.show()
fig = plt.figure()
levels = np.arange(0,.4,0.05)
plt.contour(X, Y, Tast,levels, cmap=plt.cm.hsv)
plt.colorbar()
plt.title("Contour plot of the temperature fundamental solution")
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, qast, cmap=plt.cm.hsv)
# Tweak the limits and add latex math labels.
ax.set_zlim(0, 1)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_zlabel(r'$T^*$')
plt.title("Flux fundamental solution")
plt.show()
fig = plt.figure()
q=plt.quiver(0,0,nx,ny,scale=3)
levels = np.arange(-.4,.4,0.05)
plt.contour(X, Y, qast,levels, cmap=plt.cm.hsv)
plt.colorbar()
plt.title("Contour plot of the flux fundamental solution")
plt.show()
