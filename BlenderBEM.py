import bpy
import sys
import os
import numpy as np
import bmesh

from bpy.types import Panel, Operator, PropertyGroup
from bpy.props import EnumProperty, PointerProperty, StringProperty


class MyProperties(PropertyGroup):
    
    my_enum : EnumProperty(
        name= "Enumerator / Dropdown",
        description= "Chooses the face of polygons for boundary conditions.",
        items= [('OP1', "T = 0", ""),
                ('OP2', "T = 1", ""),
                ('OP3', "q = 0", ""),
                ('OP4', "q = 1", "")
        ]
    )
    face_index = []
    face_index1 = []
    face_index2 = []
    face_index3 = []
    coord = []
    elem = []
    BCFace = []
    k = 1
    PONTOS_int = []

class ADDONNAME_PT_main_panel(Panel):
    bl_label = "BlenderBEM"
    bl_idname = "ADDONNAME_PT_main_panel"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "BlenderBEM"

    def draw(self, context):
        layout = self.layout
        scene = context.scene
        mytool = scene.my_tool
        
        layout.operator("addonname.prepare_operator")
        layout.prop(mytool, "my_enum", expand = True)
        layout.operator("addonname.submit_operator")
        layout.operator("addonname.run_operator")

class ADDONNAME_OT_prepare(Operator):
    bl_label = "Prepare mesh"
    bl_idname = "addonname.prepare_operator"
    surf = []
    elem = []
    coord = []
    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool
        enum = mytool.my_enum
        obj=bpy.context.object
        bpy.ops.object.modifier_add(type='TRIANGULATE')
        bpy.ops.object.modifier_apply(modifier="Triangulate")

        me = bpy.context.object.data
        faces = me.polygons
        verts = me.vertices
        coord = MyProperties.coord
        elem = MyProperties.elem
        for poly in me.polygons:
            elem.append([poly.vertices[0], poly.vertices[1], poly.vertices[2]])

        for v in verts:
            coord.append([v.co[0], v.co[1], v.co[2]])

        coord = np.array(coord)
        elem = np.array(elem)
        surf = np.ones(len(elem))
        ADDONNAME_OT_prepare.coord = coord
        ADDONNAME_OT_prepare.elem = elem
        ADDONNAME_OT_prepare.surf = surf
        print(f'{len(elem)} elements. Starting BEM...')
        # Visualization
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
        bpy.ops.object.mode_set(mode='OBJECT') # Can't assign materials in editmode

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


class ADDONNAME_OT_submit(Operator):
    bl_label = "Submit polys to temperature"
    bl_idname = "addonname.submit_operator"
    
    surf =[]
    elem =[]
    BCFace = []
    elemj = []
    def execute(self, context):
        scene = context.scene
        mytool = scene.my_tool
        enum = mytool.my_enum
        obj=bpy.context.object
        mats = bpy.context.object.material_slots
        face_index = mytool.face_index
        face_index1 = mytool.face_index1
        face_index2 = mytool.face_index2
        face_index3 = mytool.face_index3
        
        if obj.mode == 'EDIT':
            bm=bmesh.from_edit_mesh(obj.data)
            for p in bm.faces:
                if p.select:
                    i = p.index
                    if enum == 'OP1':
                        face_index.append(p.index)
                        mats[i].material.diffuse_color = 0, 0, 1, 1  
                    elif enum == 'OP2':
                        face_index1.append(p.index)
                        mat = obj.material_slots[i].material
                        mat.diffuse_color = 1, 0, 0, 1
                    elif enum == 'OP3':
                        face_index2.append(p.index)
                        mat = obj.material_slots[i].material
                        mat.diffuse_color = 0, 0, 0, 1
                    elif enum == 'OP4':
                        face_index3.append(p.index)
                        mat = obj.material_slots[i].material
                        mat.diffuse_color = 0, 1, 0, 1
                    
        else:
            print("Object is not in edit mode.")
        surf = ADDONNAME_OT_prepare.surf
        for i in face_index:
            surf[i] = 2

        for i in face_index1:
            surf[i] = 3
        ADDONNAME_OT_submit.surf = surf
        print(surf)
        elemj = []
        elem = ADDONNAME_OT_prepare.elem
        for i in range(len(elem)):
            elemj.append([elem[i,0],elem[i,1],elem[i,2],surf[i]])

        ADDONNAME_OT_submit.elem = elem
        ADDONNAME_OT_submit.elemj = elemj
        BCFace = []
        for i in range(1,6):
            BCFace.append([i,1,0])

        BCFace[1] = [2,0,0]
        BCFace[2] = [3,0,1]
        BCFace[3] = [4,1,0]
        BCFace[4] = [5,1,1]
        ADDONNAME_OT_submit.BCFace = BCFace

        return {'FINISHED'}


class ADDONNAME_OT_run(Operator):
    bl_label = "Run BEM"
    bl_idname = "addonname.run_operator"
    
    def execute(self, context):
        from julia.api import Julia
        jl = Julia(compiled_modules=False)
        from julia import Main

        os.chdir(os.path.dirname(bpy.data.filepath) + '/BB')
        print(os.getcwd())
        jl.eval('pwd()')
        jl.eval('include("BEM_base.jl")')

        # Solving with Julia
        # info,PONTOS_int,BCFace,k
        # NOS_GEO,ELEM,elemint,CDC = info
        Main.NOS_GEO1 = ADDONNAME_OT_prepare.coord
        Main.ELEM1 = ADDONNAME_OT_submit.elemj
        Main.elemint = []
        Main.BCFace1 = ADDONNAME_OT_submit.BCFace
        Main.CDC = []
        Main.k = 1
        Main.PONTOS_int = []
        print(Main.ELEM1)
        jl.eval('ELEM=zeros(Int,size(ELEM1,1),5)')
        jl.eval('for i=1:size(ELEM1,1); ELEM[i,:] = [i ELEM1[i,1]+1 ELEM1[i,2]+1 ELEM1[i,3]+1 ELEM1[i,4] ]; end')
        jl.eval('ELEM = convert(Array{Int32},ELEM)')
        jl.eval('NOS_GEO=zeros(size(NOS_GEO1,1),4)')
        jl.eval('for i=1:size(NOS_GEO1,1); NOS_GEO[i,:] = [i NOS_GEO1[i,1] NOS_GEO1[i,2] NOS_GEO1[i,3]]; end')

        jl.eval('BCFace=zeros(size(BCFace1,1),3)')
        jl.eval('for i=1:size(BCFace1,1); BCFace[i,:] = [BCFace1[i,1] BCFace1[i,2] BCFace1[i,3]]; end')


        T,q,T_pint,q_pint,NOS = jl.eval('potconst3d.solve([NOS_GEO,ELEM,elemint,CDC],PONTOS_int,BCFace,k)')

        T = T.real
        # Visualization
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
        bpy.ops.object.mode_set(mode='OBJECT') # Can't assign materials in editmode
        obj=bpy.context.object
        me = bpy.context.object.data
        for poly in me.polygons:
            i = poly.index
            mat = obj.material_slots[i].material
            mat.diffuse_color = 254*T[i]+1, 0, 0, 1
        
        return {'FINISHED'}


classes = [MyProperties, ADDONNAME_PT_main_panel, ADDONNAME_OT_submit,ADDONNAME_OT_prepare,ADDONNAME_OT_run]
 
def register():
    for cls in classes:
        bpy.utils.register_class(cls)
        
    bpy.types.Scene.my_tool = PointerProperty(type= MyProperties)
 
def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)
    del bpy.types.Scene.my_tool
 
if __name__ == "__main__":
    register()