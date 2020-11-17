#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 08:20:19 2017

@author: eder
"""

import boundary
import numpy as np
import graphics_problem


#extremes=np.array([-20,150])
#T=boundary.correct_T(T,extremes)

extremes=np.array([-10,3000])
graphics_problem.plot_extremes(extremes,elem,coord,T)
#
#nsurf=np.array([31])
#elsurf=graphics_problem.plot_surface(nsurf,elem,coord,surf)
#print(T[elsurf])
#
