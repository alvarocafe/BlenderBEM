#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 21:57:32 2017

@author: eder
"""

import numpy as np
def compute_shapefun(xi,eta):
    zeta=1-xi-eta
    N=np.zeros(6)
    N[0] = xi*(2*xi-1)
    N[1]=  eta*(2*eta-1)    
    N[2] = zeta*(2*zeta-1)
    N[3] = 4*xi*eta   
    N[4] = 4*eta*zeta
    N[5] = 4*xi*zeta

    return N

def compute_dshapefun(xi,eta):
    dN=np.zeros((2,6))
    zeta=1-xi-eta
    dN[0,0] = 4*xi-1
    dN[1,0] = 0
    dN[0,1] = 0
    dN[1,1] = 4*eta-1
    dN[0,2] = 1-4*zeta
    dN[1,2] = 1-4*zeta
    dN[0,3] = 4*eta
    dN[1,3] = 4*xi
    dN[0,4] = -4*eta
    dN[1,4] = 4*(zeta-eta)
    dN[0,5] = 4*(zeta-xi)
    dN[1,5] = -4*xi
    return dN