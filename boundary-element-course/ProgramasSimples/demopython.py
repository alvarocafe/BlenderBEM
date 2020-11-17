#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:40:45 2017

@author: eder
"""


def fat1(n):
    fat=1
    if(n>2):
        for i in range(1,n+1):
            fat=fat*i
    return fat

def fat2(n):
    if(n>1):
        fat=n*fat2(n-1)
    else:
        fat=1
    return fat

def fibonacci1(n):
    fib0=0
    fib1=1
    if(n==1):
        fib=0
    elif(n==2):
        fib=1
    else:
        for i in range(2,n):
            fib=fib0+fib1
            fib0=fib1
            fib1=fib
    return fib

def fibonacci2(n):
    if(n==1):
        fib=0
    elif(n==2):
        fib=1
    else:
        fib=fibonacci2(n-1)+fibonacci2(n-2)
    return fib
