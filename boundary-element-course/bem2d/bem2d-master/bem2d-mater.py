# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 20:49:28 2016
@author: mansour
Reviewed by eder on Thursday July 27 2017
"""
import numpy as np
import geometry
import solve 
import index
import input_data
import graphics_problem


bctype,bcvalue,k,file=input_data.input_data3()

inode_bound,inode_int,inode_all,coord,elem,segments,tri=geometry.compute_inodes(file)

node_med, normal = geometry.comp_node_and_normal(elem, coord)

H, G = index.mount_matrix(node_med,normal,coord,elem,k)

bcs = index.mount_bcs(segments,bctype,bcvalue)
 
A, b = index.mount_linear_system(H, G, bcs)

x = np.linalg.solve(A, b)
T, q = index.mount_vector(x, bcs)
Hin, Gin = solve.int_point(inode_int,normal,coord,elem,k)
Tint = np.dot(Hin,T) - np.dot(Gin,q)

graphics_problem.show_problem(node_med,normal,coord,bcs,tri)
graphics_problem.show_results(inode_all,inode_int,node_med,elem,coord,T,q,Tint)

# Para fazer o exercício da lista 3.
node_int_coord=np.array([[0.25,.25],[.5,.5]])
Hin2,Gin2=solve.int_point2(node_int_coord,normal,coord,elem,k)
Tint2= np.dot(Hin2,T) - np.dot(Gin2,q)
