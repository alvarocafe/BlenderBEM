#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 21:57:32 2017

@author: eder
"""

import numpy as np
def compute_shapefun(xi,eta):
    N=np.zeros(8)
    N[0] = -0.25 * (1-xi) * (1-eta) * (xi + 1+eta)
    N[4] = 0.5 * (1+xi) * (1-xi) * (1-eta)
    N[1]= 0.25 * (1+xi) * (1-eta) * (xi - 1-eta)
    N[5] = 0.5 * (1+xi) * (1+eta) * (1-eta)
    N[2] = 0.25 * (1+xi) * (1+eta) * (eta - 1+xi)
    N[6] = 0.5 * (1+eta) * (1+xi) * (1-xi)
    N[3] = 0.25 * (1-xi) * (1+eta) * (eta - 1-xi)
    N[7] = 0.5 * (1-xi) * (1+eta) * (1-eta)

    return N

def compute_dshapefun(xi,eta):
    dN=np.zeros((2,8))
    dN[0,0] = 0.25 * (1-eta) * (2 * xi + eta)
    dN[0,4] = -xi * (1-eta)
    dN[0,1] = 0.25 * (1-eta) * (2 * xi - eta)
    dN[0,5] = 0.5 * (1+eta) * (1-eta)
    dN[0,2] = 0.25 * (1+eta) * (2 * xi + eta)
    dN[0,6] = -xi * (1+eta)
    dN[0,3] = 0.25 * (1+eta) * (2 * xi - eta)
    dN[0,7] = -0.5 * (1+eta) * (1-eta)
    dN[1,0] = 0.25 * (1-xi) * (2 * eta + xi)
    dN[1,4] = -0.5 * (1+xi) * (1-xi)
    dN[1,1] = 0.25 * (1+xi) * (2 * eta - xi)
    dN[1,5] = -eta * (1+xi)
    dN[1,2] = 0.25 * (1+xi) * (2 * eta + xi)
    dN[1,6] = 0.5 * (1+xi) * (1-xi)
    dN[1,3] = 0.25 * (1-xi) * (2 * eta - xi)
    dN[1,7] = -eta * (1-xi)
    return dN