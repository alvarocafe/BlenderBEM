import bpy
import sys
import os
import numpy as np
import bmesh

from bpy.types import Panel, Operator, PropertyGroup
from bpy.props import EnumProperty, PointerProperty, StringProperty
from mathutils import Color, Vector

from julia.api import Julia
jl = Julia(compiled_modules=False)
from julia import Main

os.chdir(os.path.dirname(bpy.data.filepath) + '/BB')
print("Loading BEM solver from " + os.getcwd())
jl.eval('pwd()')
jl.eval('include("BEM_base.jl")')
print("BEM solver loaded.")

class BEMProperties(PropertyGroup):
    my_string : bpy.props.StringProperty(name= "Console")
#    my_float_vector : bpy.props.FloatVectorProperty(name= "Scale", soft_min= 0, soft_max= 1000, default= (1,1,1))
    k_float: bpy.props.FloatProperty(name = "Wavenumber (k)",default=1,min=0)
    p_float: bpy.props.FloatProperty(name = "Value",default=1)
    q_float: bpy.props.FloatProperty(name = "Flux ")
    BC_enum : EnumProperty(
        name= "Enumerator / Dropdown",
        description= "Chooses the face of polygons for boundary conditions.",
        items= [('OP1', "Potential", ""),
                ('OP2', "Flux", "")
        ]
    )
    face_index = []
    face_index1 = []

    coord = []
    elem = []
    BCFace = []
    k = 1
    PONTOS_int = []

class BLENDERBEM_PT_main_panel(Panel):
    bl_label = "BlenderBEM"
    bl_idname = "BLENDERBEM_PT_main_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "BlenderBEM"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        mytool = scene.my_tool
        
        layout.operator("blenderbem.addsphere_operator")
        layout.operator("blenderbem.prepare_operator")
        layout.prop(mytool, "k_float")
        layout.prop(mytool, "BC_enum", expand = True)
        layout.prop(mytool, "p_float")
        
        layout.operator("blenderbem.submit_operator")
        layout.operator("blenderbem.domain_operator")
        layout.operator("blenderbem.run_operator")
        layout.operator("blenderbem.runwave_operator")
#        layout.prop(mytool, "my_string")
#        layout.operator("blenderbem.julia_operator")
#        layout.operator("blenderbem.python_operator")


class BLENDERBEM_OT_JULIA(Operator):
    bl_label = "Run Julia"
    bl_idname = "blenderbem.julia_operator"

    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool
        text = mytool.my_string
        jl.eval(text)        
        return {'FINISHED'}

class BLENDERBEM_OT_PYTHON(Operator):
    bl_label = "Run Python"
    bl_idname = "blenderbem.python_operator"

    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool
        text = mytool.my_string
        eval(text)        
        return {'FINISHED'}

class BLENDERBEM_OT_prepare(Operator):
    bl_label = "Prepare mesh"
    bl_idname = "blenderbem.prepare_operator"
    surf = []
    elem = []
    coord = []
    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool
        enum = mytool.BC_enum
        obj=bpy.context.object
        bpy.ops.object.modifier_add(type='TRIANGULATE')
        bpy.ops.object.modifier_apply(modifier="Triangulate")

        me = bpy.context.object.data
        faces = me.polygons
        verts = me.vertices
        coord = BEMProperties.coord
        elem = BEMProperties.elem
        for poly in me.polygons:
            elem.append([poly.vertices[0], poly.vertices[1], poly.vertices[2]])

        for v in verts:
            coord.append([v.co[0], v.co[1], v.co[2]])

        coord = np.array(coord)
        elem = np.array(elem)
        surf = np.ones(len(elem))
        BLENDERBEM_OT_prepare.coord = coord
        BLENDERBEM_OT_prepare.elem = elem
        BLENDERBEM_OT_prepare.surf = surf
        print(f'{len(elem)} elements. Starting BEM...')
        # Visualization
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
#        bpy.ops.object.mode_set(mode='OBJECT') # Can't assign materials in editmode
        for m in me.materials:
            me.materials.pop()

        for poly in me.polygons:
            i = poly.index
            # TODO: check if the object already has materials.
            mat = bpy.data.materials.new("Mat_%i" % i)
        #    mat.use_diffuse_ramp = True
        #    mat.diffuse_ramp.evaluate(T)
            mat.diffuse_color = 255, 255, 255, 1
            mat.roughness = 1
            me.materials.append(mat)
            poly.material_index = i

        return {'FINISHED'}


class BLENDERBEM_OT_submit(Operator):
    bl_label = "Submit polys to boundary condition"
    bl_idname = "blenderbem.submit_operator"
    surf = []
    elem =[]
    elemj = []
    BCFace = [[1,1,0.0],[2,0,0.0]]
    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool
        enum = mytool.BC_enum
        pvalue = mytool.p_float        
        obj=bpy.context.object
        mats = bpy.context.object.material_slots
        BCFace = BLENDERBEM_OT_submit.BCFace
        surf = BLENDERBEM_OT_prepare.surf
        if obj.mode == 'EDIT':
            bm=bmesh.from_edit_mesh(obj.data)
            for p in bm.faces:
                if p.select:
                    i = p.index
                    if enum == 'OP1':
                        for BC in BCFace:
                            if pvalue == BC[2]:
                                surf[i] = BC[0]
#                            elif enum == 'OP1':
                            else:
                                surf[p.index] = max(surf) + 1
                                BCFace.append([max(surf),0,pvalue])
                        mats[i].material.diffuse_color = 0, 0, 1, 1
                    elif enum == 'OP2':
                        for BC in BCFace:
                            if pvalue == BC[2]:
                                surf[i] = BC[0]
#                            elif enum == 'OP1':
                            else:
                                surf[p.index] = max(surf) + 1
                                BCFace.append([max(surf),1,pvalue])
                        mats[i].material.diffuse_color = 1, 0, 0, 1
        else:
            print("Object is not in edit mode.")

        BLENDERBEM_OT_submit.surf = surf
        BLENDERBEM_OT_prepare.BCFace = BCFace
        
        elemj = []
        elem = BLENDERBEM_OT_prepare.elem
        for i in range(len(elem)):
            elemj.append([elem[i,0],elem[i,1],elem[i,2],surf[i]])
        BLENDERBEM_OT_submit.elem = elem
        BLENDERBEM_OT_submit.elemj = elemj
        BLENDERBEM_OT_submit.BCFace = BCFace

        return {'FINISHED'}

class BLENDERBEM_OT_domain(Operator):
    bl_label = "Submit object to domain"
    bl_idname = "blenderbem.domain_operator"
    object = ''
    PONTOS_dom = []
    def execute(self, context):
        
        PONTOS_dom = []
        me = bpy.context.object.data
        bpy.ops.object.modifier_add(type='TRIANGULATE')
        bpy.ops.object.modifier_apply(modifier="Triangulate")        
        BLENDERBEM_OT_domain.object = bpy.context.object
        faces = me.polygons
        verts = me.vertices
        coord = []
        elem = []
        ndom = 0
        for poly in faces:
            elem.append([poly.index,poly.vertices[0], poly.vertices[1], poly.vertices[2],1])

        for v in verts:
            coord.append([v.index+1,v.co[0], v.co[1], v.co[2]])
        Main.NOS_GEO_dom = coord
        Main.ELEM_dom = elem
        
#        PONTOS_dom = jl.eval('const3D_tri.mostra_geoTRI(NOS_GEO_dom,ELEM_dom)')
        PONTOS_dom = coord

        BLENDERBEM_OT_domain.PONTOS_dom = PONTOS_dom
        print(f'{len(PONTOS_dom)} domain points...')
        # Visualization
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
#        for poly in me.polygons:
#            i = poly.index
#            # TODO: check if the object already has materials.
#            mat = bpy.data.materials.new("Mat_%i" % i)
#        #    mat.use_diffuse_ramp = True
#        #    mat.diffuse_ramp.evaluate(T)
#            mat.diffuse_color = 255, 255, 255, 1
#            mat.roughness = 1
#            me.materials.append(mat)
#            poly.material_index = i

        return {'FINISHED'}



class BLENDERBEM_OT_run(Operator):
    bl_label = "Run Laplace BEM"
    bl_idname = "blenderbem.run_operator"

    def color_verts(my_object,T):
        vert_list = my_object.vertices
        color_map = my_object.vertex_colors.new()
        i = 0
        for poly in my_object.polygons:
            for idx in poly.loop_indices:
                loop = my_object.loops[idx]
                v = loop.vertex_index
                color_map.data[i].color = 0, 0, 255*T[v], 1
                i += 1
        return {'FINISHED'}

    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool
        # Solving with Julia
        # info,PONTOS_int,BCFace,k
        # NOS_GEO,ELEM,elemint,CDC = info
        Main.NOS_GEO1 = BLENDERBEM_OT_prepare.coord
        Main.ELEM1 = BLENDERBEM_OT_submit.elemj
        Main.elemint = []
        Main.BCFace1 = BLENDERBEM_OT_submit.BCFace
        Main.CDC = []
        Main.k = mytool.k_float
        Main.PONTOS_int = BLENDERBEM_OT_domain.PONTOS_dom
        PONTOS_int = BLENDERBEM_OT_domain.PONTOS_dom
        object = BLENDERBEM_OT_domain.object
        jl.eval('ELEM=zeros(Int,size(ELEM1,1),5)')
        jl.eval('for i=1:size(ELEM1,1); ELEM[i,:] = [i ELEM1[i,1]+1 ELEM1[i,2]+1 ELEM1[i,3]+1 ELEM1[i,4] ]; end')
        jl.eval('ELEM = convert(Array{Int32},ELEM)')
        jl.eval('NOS_GEO=zeros(size(NOS_GEO1,1),4)')
        jl.eval('for i=1:size(NOS_GEO1,1); NOS_GEO[i,:] = [i NOS_GEO1[i,1] NOS_GEO1[i,2] NOS_GEO1[i,3]]; end')

        jl.eval('BCFace=zeros(size(BCFace1,1),3)')
        jl.eval('for i=1:size(BCFace1,1); BCFace[i,:] = [BCFace1[i,1] BCFace1[i,2] BCFace1[i,3]]; end')


        T,q,T_pint,q_pint,NOS = jl.eval('potconst3d.solve([NOS_GEO,ELEM,elemint,CDC],PONTOS_int,BCFace,k)')

        T = T.real
        T_pint = T_pint.real

        # Visualization
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
        obj=bpy.context.object
        me = bpy.context.object.data
        for poly in me.polygons:
            i = poly.index
            mat = obj.material_slots[i].material
            mat.diffuse_color = 254*T[i]+1, 0, 0, 1

        # Color domain points with vertex color
        if BLENDERBEM_OT_domain.object != '':
#            print(object)
            BLENDERBEM_OT_run.color_verts(object,T_pint)
#            object = bpy.data.objects[object]           
        return {'FINISHED'}



class BLENDERBEM_OT_runwave(Operator):
    bl_label = "Run Helmholtz BEM"
    bl_idname = "blenderbem.runwave_operator"

    def color_verts(my_object,T):
        my_object = my_object.data
        vert_list = my_object.vertices
        color_map = my_object.vertex_colors.new()
        i = 0
        for poly in my_object.polygons:
            for idx in poly.loop_indices:
                loop = my_object.loops[idx]
                v = loop.vertex_index
                color_map.data[i].color = 0, 0, 254*T[v]+1, 1
                i += 1
        return {'FINISHED'}
    
    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool        
        # Solving with Julia
        # info,PONTOS_int,BCFace,k
        # NOS_GEO,ELEM,elemint,CDC = info
        Main.NOS_GEO1 = BLENDERBEM_OT_prepare.coord
        Main.ELEM1 = BLENDERBEM_OT_submit.elemj
        Main.elemint = []
        Main.BCFace1 = BLENDERBEM_OT_submit.BCFace
        Main.CDC = []
        Main.k = mytool.k_float
        Main.PONTOS_int = BLENDERBEM_OT_domain.PONTOS_dom
        PONTOS_int = BLENDERBEM_OT_domain.PONTOS_dom
        object = BLENDERBEM_OT_domain.object
        jl.eval('ELEM=zeros(Int,size(ELEM1,1),5)')
        jl.eval('for i=1:size(ELEM1,1); ELEM[i,:] = [i ELEM1[i,1]+1 ELEM1[i,2]+1 ELEM1[i,3]+1 ELEM1[i,4] ]; end')
        jl.eval('ELEM = convert(Array{Int32},ELEM)')
        jl.eval('NOS_GEO=zeros(size(NOS_GEO1,1),4)')
        jl.eval('for i=1:size(NOS_GEO1,1); NOS_GEO[i,:] = [i NOS_GEO1[i,1] NOS_GEO1[i,2] NOS_GEO1[i,3]]; end')

        jl.eval('BCFace=zeros(size(BCFace1,1),3)')
        jl.eval('for i=1:size(BCFace1,1); BCFace[i,:] = [BCFace1[i,1] BCFace1[i,2] BCFace1[i,3]]; end')


        T,q,T_pint,q_pint = jl.eval('const3D_tri.solve([NOS_GEO,ELEM,elemint,CDC],PONTOS_int,BCFace,k,true)')

        T = T.real
        T_pint = T_pint.real
#        T_pint = abs(T_pint/max(T_pint))
        # Visualization
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
        obj=bpy.context.object
        me = bpy.context.object.data
        for poly in me.polygons:
            i = poly.index
            mat = obj.material_slots[i].material
            mat.diffuse_color = 0, 0, 254*T[i]+1, 1
        # Color domain points with materials
        if BLENDERBEM_OT_domain.object != '':
            obj=BLENDERBEM_OT_domain.object
            me = BLENDERBEM_OT_domain.object.data
#            for poly in me.polygons:
#                i = poly.index
#                mat = obj.material_slots[i].material
#                mat.diffuse_color = 0, 0, 254*T_pint[i]+1, 1            
            BLENDERBEM_OT_runwave.color_verts(obj,T_pint)
        return {'FINISHED'}

class BLENDERBEM_OT_addcube(Operator):
    bl_label = "Add an exterior cube"
    bl_idname = "blenderbem.addcube_operator"

    def execute(self, context):
        bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        bpy.ops.object.subdivision_set(level=3, relative=False)
        bpy.context.object.modifiers["Subdivision"].subdivision_type = 'SIMPLE'
        bpy.ops.object.modifier_apply(modifier="Subdivision")
        bpy.ops.object.editmode_toggle()
        bpy.ops.mesh.normals_make_consistent(inside=True)
        bpy.ops.object.editmode_toggle()

        return {'FINISHED'}
class BLENDERBEM_OT_addsphere(Operator):
    bl_label = "Add an exterior sphere"
    bl_idname = "blenderbem.addsphere_operator"

    def execute(self, context):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=1, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        bpy.ops.object.editmode_toggle()
        bpy.ops.mesh.normals_make_consistent(inside=True)
        bpy.ops.object.editmode_toggle()
        bpy.ops.mesh.primitive_plane_add(size=2, enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        bpy.ops.transform.resize(value=(10, 10, 10), orient_type='GLOBAL', orient_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)), orient_matrix_type='GLOBAL', mirror=True, use_proportional_edit=False, proportional_edit_falloff='SMOOTH', proportional_size=1, use_proportional_connected=False, use_proportional_projected=False)
        bpy.ops.object.modifier_add(type='SUBSURF')
        bpy.context.object.modifiers["Subdivision"].subdivision_type = 'SIMPLE'
        bpy.context.object.modifiers["Subdivision"].levels = 6
        bpy.ops.object.modifier_apply(modifier="Subdivision")

        return {'FINISHED'}


classes = [BEMProperties, BLENDERBEM_PT_main_panel, BLENDERBEM_OT_submit,BLENDERBEM_OT_prepare,BLENDERBEM_OT_run,BLENDERBEM_OT_runwave,BLENDERBEM_OT_domain,BLENDERBEM_OT_addcube,BLENDERBEM_OT_addsphere, BLENDERBEM_OT_JULIA, BLENDERBEM_OT_PYTHON]
 
def register():
    for cls in classes:
        bpy.utils.register_class(cls)
        
    bpy.types.Scene.my_tool = PointerProperty(type= BEMProperties)
 
def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)
    del bpy.types.Scene.my_tool
 
if __name__ == "__main__":
    register()