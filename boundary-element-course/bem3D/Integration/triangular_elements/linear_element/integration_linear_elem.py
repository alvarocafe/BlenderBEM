# -*- coding: utf-8 -*-
"""
Spyder Editor

Este é um arquivo de script temporário.
"""
import numpy as np
import shape_functions as sp
import matplotlib.pyplot as plt
plt.close("all")
x=np.zeros(9)
y=np.zeros(9)
x[0]=0
y[0]=0
x[1]=2
y[1]=0
x[2]=2
y[2]=3
x[3]=0
y[3]=2.

elem=np.array([[0,1,2],[0,2,3]],dtype=int)
nelem=elem.shape[0]

plt.figure()
ngp=4
xiquad,wquad=np.polynomial.legendre.leggauss(ngp)
a=0
b=1
xitri = 0.5*(xiquad + 1)*(b - a) + a
wtri=wquad*0.5*(b - a)

eta=xitri
I=0
for k in range(nelem):
    nos=elem[k,:]
    xnodes=x[nos]
    ynodes=y[nos]
    for j in range(ngp):
        for i in range(ngp):
            xi=(1-eta[j])*xitri[i]
            N=sp.compute_shapefun(xi,eta[j])
            dN=sp.compute_dshapefun(xi,eta[j])
            xx=np.dot(N,xnodes)
            yy=np.dot(N,ynodes)
            dxdxi=np.dot(dN[0,:],xnodes)
            dxdeta=np.dot(dN[1,:],xnodes)
            dydxi=np.dot(dN[0,:],ynodes)
            dydeta=np.dot(dN[1,:],ynodes)
            J=dxdxi*dydeta-dxdeta*dydxi
            f=(xx**2+yy);
            I=I+f*(1-eta[j])*J*wtri[i]*wtri[j]
            plt.plot(np.append(xnodes[0:3],xnodes[0]),np.append(ynodes[0:3],ynodes[0]),'ro-')
            plt.plot(xnodes,ynodes,'gd')
print('The value of the integral is: ',I)
ax = plt.axes()
for i in range(4):
    ax.text(x[i], y[i], i+1, fontsize=15)
plt.show()
