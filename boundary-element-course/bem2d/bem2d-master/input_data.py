#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:08:57 2017

@author: eder
"""
def input_data1():
    # boundary type [segm1 segm1 ...]
    # types: 0: dirichlet; 1: neumann
    bt = [1, 0, 1, 0]
    # boundary condition [segm1 segm2 ...]
    bc = [2, 1 , 0, 0]
    k=1
    file ='squaremesh0'
    return bt,bc,k,file

def input_data2():
    # boundary type [segm1 segm1 ...]
    # types: 0: dirichlet; 1: neumann
    btype = [0, 0, 1, 0, 0, 1, 1]
    # boundary condition [segm1 segm2 ...]
    bc = [1, 2, 3, 0, 5, 1, 1]
    k=1
    file = 'placa_furo'
    return btype,bc,k,file

def input_data3():
    # boundary type [segm1 segm1 ...]
    # types: 0: dirichlet; 1: neumann
    btype = [0, 0, 1, 0, 0, 1, 1,1,1]
    # boundary condition [segm1 segm2 ...]
    bc = [1, 2, 3,0, 5, 1, 1,1,1]
    k=1
    file = 'placa_furo2'
    return btype,bc,k,file

