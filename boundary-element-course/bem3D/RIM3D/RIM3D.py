# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 22:29:08 2016
example of 3d bem program using the c++ solver
@author: mansour
require normals pointing outward.
Reviewed by Eder on Sunday, July 30rd, 18:07 2017
"""
from mesh import Mesh
import RIM_functions
import numpy as np
import graphics_problem
import input_data

file,surf_bc,k=input_data.input_data5()

mymesh = Mesh() #criou objeto: instanciou a classe
#chamei o metodo do objeto
mymesh.read_msh(file + '.msh')
coord = mymesh.Verts
# get only the triangular elements [[physid],[nodes]] physid is the surf of the elem

elem = mymesh.Elmts[2][1]-1
surf = mymesh.Elmts[2][0]-1

# calculate the mean node and normal
normal, Jac = RIM_functions.compute_normal(coord, elem)

surface,volume = RIM_functions.compute_geoprop(normal,Jac,coord, elem,k)

graphics_problem.show_results(elem,coord)
print('Surface area = ',surface)
print('Volume = ',volume)
print('Number of elements = ',elem.shape[0])
