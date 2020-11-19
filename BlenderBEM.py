import bpy
import sys
import os

dir = os.path.dirname(bpy.data.filepath)
dir += '/boundary-element-course/bem3D/bem3d-master/'
if not dir in sys.path:
    sys.path.append(dir)

os.chdir(os.path.dirname(bpy.data.filepath) + '/boundary-element-course/bem3D/bem3d-master/')
print(os.getcwd())

from mesh import Mesh
import boundary
import geometry
import index
import numpy as np
import graphics_problem
import input_data

bpy.ops.mesh.primitive_cube_add(size=2000, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
bpy.ops.object.modifier_add(type='TRIANGULATE')
bpy.ops.object.modifier_apply(modifier="Triangulate")

me = bpy.context.object.data
uv_layer = me.uv_layers.active.data
faces = me.polygons
verts = me.vertices
coord = []
elem = []

for poly in me.polygons:
    print("Polygon index: %d, length: %d" % (poly.index, poly.loop_total))
    print("Vertices in polygon %d:" % poly.index)
    print(list(poly.vertices))
    elem.append([poly.vertices[0], poly.vertices[1], poly.vertices[2]])

for v in verts:
    print("Vertex: %d" % v.index)
    coord.append([v.co[0], v.co[1], v.co[2]])

coord = np.array(coord)
elem = np.array(elem)
surf = np.array([0,0,1,1,2,2,3,3,4,4])
surf_bc={0:[0,0],1:[0,1]}
k = 1
coord_med, normal, Jac = geometry.node_med(coord, elem)

H, G = index.mount_matrix(coord_med,normal,Jac,coord, elem,k)

# apply boundary condition
elem_bc = boundary.elem(surf_bc, surf, elem)
 
A, b = index.mount_linear_system(H, G, elem_bc)

x = np.linalg.solve(A, b)
T, q = index.mount_vector(x, elem_bc)