# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 22:29:08 2016
example of 3d bem program using the c++ solver
@author: mansour
require normals pointing outward.
Reviewed by Eder on Sunday, July 30rd, 18:07 2017
"""
from mesh import Mesh
import boundary
import geometry
import index
import numpy as np
import graphics_problem
import input_data


file,surf_bc,k=input_data.input_data4()

mymesh = Mesh() #criou objeto: instanciou a classe
#chamei o metodo do objeto
mymesh.read_msh(file + '.msh')
coord = mymesh.Verts
# get only the triangular elements [[physid],[nodes]] physid is the surf of the elem

elem = mymesh.Elmts[2][1]-1
surf = mymesh.Elmts[2][0]-1



# calculate the mean node and normal
coord_med, normal, Jac = geometry.node_med(coord, elem)

print(coord_med)

H, G = index.mount_matrix(coord_med,normal,Jac,coord, elem,k)

# apply boundary condition
elem_bc = boundary.elem(surf_bc, surf, elem)

A, b = index.mount_linear_system(H, G, elem_bc)

print(H)
print(G)

x = np.linalg.solve(A, b)
T, q = index.mount_vector(x, elem_bc)
graphics_problem.show_results(elem,coord,T)
