#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 10:36:07 2017

@author: eder
"""

import numpy as np
import matplotlib.pyplot as plt

def compute_dshapefun(xi,eltype):
    if(eltype==1):
        dN1dxi=-0.5 + xi
        dN2dxi=-2*xi
        dN3dxi=0.5 + 1.*xi
    else:
        dN1dxi=3*(-1 + 3*xi)/4
        dN2dxi=-9*xi/2
        dN3dxi=3*(1 + 3*xi)/4
    return  dN1dxi,dN2dxi,dN3dxi

def compute_shapefun(xi,eltype):
    if(eltype==1):
        N1 = 0.5*xi*(xi-1)
        N2 = 1-xi**2
        N3 = 0.5*xi*(xi+1)
    else:
        N1 = xi*(9/8*xi-3/4)
        N2 = 1-9/4*xi**2
        N3 = xi*(9/8*xi+3/4)
    return N1,N2,N3



    
R=1

eltype=2; # eltype = 1 = > continuous, eltype = 2 => discontinuous

if(eltype==2):
    thetas=np.array([-15,-45,-75])
else:
    thetas=np.array([0,-45,-90])
    


x1=R*np.cos(thetas[0]*np.pi/180)
x2=R*np.cos(thetas[1]*np.pi/180)
x3=R*np.cos(thetas[2]*np.pi/180)
y1=R*np.sin(thetas[0]*np.pi/180)
y2=R*np.sin(thetas[1]*np.pi/180)
y3=R*np.sin(thetas[2]*np.pi/180)
L=np.pi/2*R # Length of the segment
vetx=np.array([x1,x2,x3]) # x node coordinates
vety=np.array([y1,y2,y3]) # y node coordinates
xi=np.linspace(-1,1,10)
ngp=len(xi)
L2=L/2*np.ones(ngp)
dgamadxi=np.zeros(ngp)
x=np.zeros(ngp)
y=np.zeros(ngp)
xcirc=np.zeros(ngp)
ycirc=np.zeros(ngp)

for i in range(ngp):
    dN1dxi,dN2dxi,dN3dxi=compute_dshapefun(xi[i],eltype)
    N1,N2,N3=compute_shapefun(xi[i],eltype)
    x[i]=N1*x1+N2*x2+N3*x3
    y[i]=N1*y1+N2*y2+N3*y3
    dxdxi=dN1dxi*x1+dN2dxi*x2+dN3dxi*x3
    dydxi=dN1dxi*y1+dN2dxi*y2+dN3dxi*y3
    dgamadxi[i]=np.sqrt(dxdxi**2+dydxi**2)
    theta=(xi[i]-1)*45
    xcirc[i]=R*np.cos(theta*np.pi/180)
    ycirc[i]=R*np.sin(theta*np.pi/180)

plt.figure()
plt.plot(x,y,'r-',label='boundary element');
plt.plot(xcirc,ycirc,'k--',label='circunference');
plt.plot(vetx,vety,'bo',label='nodes');
plt.axis("equal")
plt.legend(loc="upper left",fancybox="true")

plt.show()
plt.figure()
plt.plot(xi,dgamadxi,'r-',label='jacobian')
plt.plot(xi,L2,'b--',label='L/2')
plt.legend(loc="upper center",fancybox="true")
plt.show()

