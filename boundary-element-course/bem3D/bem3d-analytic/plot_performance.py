#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 10:06:06 2017

@author: eder
"""

import numpy as np
import matplotlib.pyplot as plt
x = np.array([612,2480,3032,9772,9792,21056])
y1 = np.array([1.7842756419995567,28.281876537000016,48.21456806199967,497.5592855200001,449.90157704599915,2058.8115712870003])
y2 = np.array([0.1293355289963074,0.3038316990000567,0.43758968400106824,15.369493176000105,17.567136429999664,287.22895112399965 ])
plt.close("all")
plt.figure()
plt.plot(x,y1,linestyle='-',marker='x')
plt.axis("tight") # Fit the axis tightly to the plot
plt.title("Time to assembly matrix A and vector b")
plt.xlabel("Nuber of nodes")
plt.ylabel("Time (seconds)")
plt.grid("on")

plt.figure()
plt.plot(x,y2,linestyle='-',marker='x')
plt.axis("tight") # Fit the axis tightly to the plot
plt.title("Time to solve the linear system A x = b")
plt.xlabel("Nuber of nodes")
plt.ylabel("Time (seconds)")
plt.grid("on")


plt.show()
