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
import cindex
import cindex2
import index
import numpy as np
import graphics_problem
import input_data
import timeit

totaltime=0
tic = timeit.default_timer()
file,surf_bc,k=input_data.input_data5()
mymesh = Mesh() #criou objeto: instanciou a classe
#chamei o metodo do objeto
mymesh.read_msh(file + '.msh')
coord = mymesh.Verts
# get only the triangular elements [[physid],[nodes]] physid is the surf of the elem
elem = mymesh.Elmts[2][1]-1
surf = mymesh.Elmts[2][0]-1
# calculate the mean node and normal
coord_med, normal, Jac = geometry.node_med(coord, elem)
# apply boundary condition
elem_bc = boundary.elem(surf_bc, surf, elem)
toc = timeit.default_timer()
print("Time to read msh file: ", toc-tic," seconds")
totaltime=totaltime+toc-tic

# In order to compile pyx files, in the terminal use the command:
# python setup.py build_ext -i

tic = timeit.default_timer()
A,b = cindex2.mount_matrix(coord_med,normal,Jac,coord, elem,k,elem_bc)
toc = timeit.default_timer()
print("Cython C time to assembly matrix A and vector b: ", toc-tic," seconds")
totaltime=totaltime+toc-tic
#
#tic = timeit.default_timer()
#A,b = cindex.mount_matrix(coord_med,normal,Jac,coord, elem,k,elem_bc)
#toc = timeit.default_timer()
#print("Cython Numpy time to assembly matrix A and vector b: ", toc-tic," seconds")
#
#tic = timeit.default_timer()
#A,b = index.mount_matrix(coord_med,normal,Jac,coord, elem,k,elem_bc)
#toc = timeit.default_timer()
#print("Python time to assembly matrix A and vector b: ", toc-tic," seconds")

tic = timeit.default_timer()
x = np.linalg.solve(A, b)
toc = timeit.default_timer()
print("Time to solve the linear system: ", toc-tic," seconds")
totaltime=totaltime+toc-tic

tic = timeit.default_timer()
T, q = index.mount_vector(x, elem_bc)
toc = timeit.default_timer()
print("Time to separate temperature and flux: ", toc-tic," seconds")
totaltime=totaltime+toc-tic


tic = timeit.default_timer()
graphics_problem.show_results(elem,coord,T)
toc = timeit.default_timer()
print("Time to post-processing: ", toc-tic," seconds")
totaltime=totaltime+toc-tic

print("Total time: ", totaltime," seconds")
print("Number of elements: ", elem.shape[0])

