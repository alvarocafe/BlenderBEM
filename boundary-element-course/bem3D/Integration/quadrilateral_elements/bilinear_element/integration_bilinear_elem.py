# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""
import numpy as np
import shape_functions as sp
import matplotlib.pyplot as plt
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
ngp=4
xi,w=np.polynomial.legendre.leggauss(ngp)
eta=xi
ro=w
I=0
for j in range(ngp):
    for i in range(ngp):
        N=sp.compute_shapefun(xi[i],eta[j])
        dN=sp.compute_dshapefun(xi[i],eta[j])
        xx=np.dot(N,x)
        yy=np.dot(N,y)
        dxdxi=np.dot(dN[0,:],x)
        dxdeta=np.dot(dN[1,:],x)
        dydxi=np.dot(dN[0,:],y)
        dydeta=np.dot(dN[1,:],y)
        J=dxdxi*dydeta-dxdeta*dydxi;
        f=(xx**2+yy);
        I=I+f*J*w[i]*ro[j]
        
print('The value of the integral is: ',I)
plt.figure()
plt.plot(np.append(x,x[0]),np.append(y,y[0]),'ro-')
ax = plt.axes()
for i in range(4):
    ax.text(x[i], y[i], i+1, fontsize=15)
plt.show()
