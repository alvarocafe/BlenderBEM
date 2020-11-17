#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:08:57 2017

@author: eder
"""
def input_data1():
    #name of the file
    #file = 'PlacaFuro';
    file = 'viga';
    #boundary condition 0:temp 1:flux, not define is zero flux [1, 0] where [type, value]
    surf_bc={0:[0,0],1:[0,1]}
    k=1
    return file,surf_bc,k

def input_data2():
    #name of the file
    #file = 'PlacaFuro';
    file = 'Placa_furo';
    #boundary condition 0:temp 1:flux, not define is zero flux [1, 0] where [type, value]
    surf_bc={3:[0,0],1:[0,1]}
    k=1
    return file,surf_bc,k

def input_data3():
    #name of the file
    #file = 'PlacaFuro';
    file = 'L';
    #boundary condition 0:temp 1:flux, not define is zero flux [1, 0] where [type, value]
    surf_bc={3:[0,0],1:[0,1]}
    k=1
    return file,surf_bc,k


def input_data4():
    #name of the file
    #file = 'PlacaFuro';
    file = 'Tutorial_v4';
    #boundary condition 0:temp 1:flux, not define is zero flux [1, 0] where [type, value]
    surf_bc={3:[0,0],1:[0,1]}
    k=1
    return file,surf_bc,k


def input_data5():
    #name of the file
    #file = 'PlacaFuro';
    file = 'Polia';
    #boundary condition 0:temp 1:flux, not define is zero flux [1, 0] where [type, value]
    surf_bc={3:[0,0],1:[0,1]}
    k=1
    return file,surf_bc,k

