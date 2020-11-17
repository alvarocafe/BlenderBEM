#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:40:58 2017

@author: eder
"""

import matplotlib.pyplot as plt
import numpy as np
plt.close("all")
x = np.linspace(0,2*np.pi,50)
y1 = np.sin(x)
y2 = np.cos(x)
p1 = plt.plot(x,y1,linestyle='-',marker='x',label='sin(x)')
p2 = plt.plot(x,y2,linestyle="--",marker="o",label="cos(x)")
plt.axis("tight") # Fit the axis tightly to the plot
plt.title("funções trigonométricas")
plt.grid("on")
# Create a legend of all the existing plots using their labels as names
plt.legend(loc="upper center",fancybox="true")
plt.show()
