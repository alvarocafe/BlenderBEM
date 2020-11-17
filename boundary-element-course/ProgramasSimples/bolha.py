#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:20:03 2017

@author: eder
"""
import numpy as np
n=10;
a=np.random.randint(1,101,size=n)
print(a)
for i in range(n-1):
    for j in range(n-i-1):
        print("i = ",i," j = ",j)
        if(a[j]>a[j+1]):
            a[j],a[j+1]=a[j+1],a[j]
		
print(a)
		  