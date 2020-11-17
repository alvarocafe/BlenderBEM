#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 15:17:11 2017

@author: eder
"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def show_results(tri,XYZ,T):
    plt.close("all")
    x=np.zeros(3)
    y=np.zeros(3)
    z=np.zeros(3)
    zc=np.zeros(len(tri))
    pc=[]
    for elem in range(len(tri)):
        no1=tri[elem,0]
        no2=tri[elem,1]
        no3=tri[elem,2]
        x[0]=XYZ[no1,0]
        y[0]=XYZ[no1,1]
        z[0]=XYZ[no1,2]
        x[1]=XYZ[no2,0]
        y[1]=XYZ[no2,1]
        z[1]=XYZ[no2,2]
        x[2]=XYZ[no3,0]
        y[2]=XYZ[no3,1]
        z[2]=XYZ[no3,2]
        pc.append([[x[0],y[0],z[0]],[x[1],y[1],z[1]],[x[2],y[2],z[2]]])
        zc[elem]=T[elem]
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    m = cm.ScalarMappable(cmap=cm.jet)
    b=m.to_rgba(zc)
    vetcor=[(i[0],i[1],i[2]) for i in b]
    m.set_array([min(zc),max(zc)])
    m.set_clim(vmin=min(zc),vmax=max(zc))
    q = Poly3DCollection(pc, linewidths=1,edgecolors="k")
    q.set_facecolor(vetcor)
    ax.add_collection3d(q)
    fig.colorbar(m)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    xmin=min(XYZ[:,0])
    ymin=min(XYZ[:,1])
    zmin=min(XYZ[:,2])
    xmax=max(XYZ[:,0])
    ymax=max(XYZ[:,1])
    zmax=max(XYZ[:,2])
    deltax=xmax-xmin
    deltay=ymax-ymin    
    deltaz=zmax-zmin
    deltamax=np.max([deltax,deltay,deltaz])
    ax.set_xlim3d(xmin-.2*deltamax,xmin+1.*deltamax)
    ax.set_ylim3d(ymin-.2*deltamax,ymin+1.*deltamax)
    ax.set_zlim3d(zmin-.2*deltamax,zmin+1.*deltamax)
    ax.view_init(elev=18., azim=43.)
    plt.title('Temperature')
    ax.axis('off')    
    plt.show()
