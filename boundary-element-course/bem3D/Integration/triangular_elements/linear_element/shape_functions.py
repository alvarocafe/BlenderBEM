#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 21:57:32 2017

@author: eder
"""

import numpy as np
def compute_shapefun(xi,eta):
    zeta=1-xi-eta
    N=np.zeros(3)
    N[0] = xi
    N[1] = eta
    N[2] = zeta

    return N

def compute_dshapefun(xi,eta):
    dN=np.zeros((2,3))
    zeta=1-xi-eta
    dN[0,0] = 1
    dN[1,0] = 0
    dN[0,1] = 0
    dN[1,1] = 1
    dN[0,2] = -1
    dN[1,2] = -1
    return dN