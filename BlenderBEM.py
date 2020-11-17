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
coords = []
elems = []
for poly in me.polygons:
    print("Polygon index: %d, length: %d" % (poly.index, poly.loop_total))
    print("Vertices in polygon %d:" % poly.index)
    print(list(poly.vertices))
    for loop_index in range(poly.loop_start, poly.loop_start + poly.loop_total):
        print("    Vertex: %d" % me.loops[loop_index].vertex_index)
        print("    UV: %r" % uv_layer[loop_index].uv)
        print("coords: ")
        print(verts[me.loops[loop_index].vertex_index].co)
        
coords = []
elems = []
for poly in me.polygons:
    print("Polygon index: %d, length: %d" % (poly.index, poly.loop_total))
    print("Vertices in polygon %d:" % poly.index)
    print(list(poly.vertices))
    elems.append([poly.vertices[0], poly.vertices[1], poly.vertices[2]])
    
for v in verts:
    print("Vertex: %d" % v.index)
    coords.append([v.co[0], v.co[1], v.co[2]])
coords = np.array(coords)
elems = np.array(elems)