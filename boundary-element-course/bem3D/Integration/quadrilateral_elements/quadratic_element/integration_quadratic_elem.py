# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""
import numpy as np
import shape_functions as sp
import matplotlib.pyplot as plt
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
ngp=4
xi,weight=np.polynomial.legendre.leggauss(ngp)
eta=xi
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
        I=I+f*J*weight[i]*weight[j]
        
print('The value of the integral is: ',I)