#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 21:57:32 2017

@author: eder
"""

import numpy as np
def compute_shapefun(xi,eta):
    N=np.zeros(4)
    N[0]=1/4*(1-xi)*(1-eta)
    N[1]=1/4*(1+xi)*(1-eta)
    N[2]=1/4*(1+xi)*(1+eta)
    N[3]=1/4*(1-xi)*(1+eta)
    return N

def compute_dshapefun(xi,eta):
    dN=np.zeros((2,4))
    dN[0,0]=-1/4*(1-eta)
    dN[1,0]=-1/4*(1-xi)
    dN[0,1]=1/4*(1-eta)
    dN[1,1]=-1/4*(1+xi)
    dN[0,2]=1/4*(1+eta)
    dN[1,2]=1/4*(1+xi)
    dN[0,3]=-1/4*(1+eta)
    dN[1,3]=1/4*(1-xi)
    return dN