# -*- coding: utf-8 -*-
"""
Universidade de Brasília
Departamento de Engenharia Mecânica
Brasília, dezembro de 2020

Programa de elementos de contorno aplicado a problemas de condução de
calor tri-dimensional sem fontes de calor concentradas

Tipo de elementos: Triangulares lineares contínuos

Autores: Éder Lima de Albuquerque
         Gustavo Silva Vaz Gontijo (ggontijo@gmail.com)
         Inácio Miura
         Álvaro Campos Ferreira (alvaro.campos.ferreira@gmail.com)

"""
#%% BIBLIOTECAS E ARQUIVOS NECESSÁRIOS
import bpy
import sys
import os
import numpy as np
import bmesh
import time

from bpy.types import Panel, Operator, PropertyGroup
from bpy.props import EnumProperty, PointerProperty, StringProperty
from mathutils import Color, Vector

dir = os.path.dirname(bpy.data.filepath)
dir += '/boundary-element-course/bem3D/potlinear3D/'
if not dir in sys.path:
    sys.path.append(dir)

os.chdir(os.path.dirname(bpy.data.filepath) + '/boundary-element-course/bem3D/potlinear3D/')
print(os.getcwd())

import entrada_de_dados
import contorno
import sistema
import integracao as integ


class MyProperties(PropertyGroup):
    
    my_enum : EnumProperty(
        name= "Enumerator / Dropdown",
        description= "Chooses the face of polygons for boundary conditions.",
        items= [('OP1', "p = 0", ""),
                ('OP2', "p = 1", ""),
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

        layout.operator("addonname.addcube_operator")        
        layout.operator("addonname.addsphere_operator")
        layout.operator("addonname.prepare_operator")
        layout.prop(mytool, "my_enum", expand = True)
        layout.operator("addonname.submit_operator")
#        layout.operator("addonname.domain_operator")
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
        surf = np.repeat('flux0',len(elem))
        ADDONNAME_OT_prepare.coord = coord
        ADDONNAME_OT_prepare.elem = elem
        ADDONNAME_OT_prepare.surf = surf
        print(f'{len(elem)} elements. Starting BEM...')
        # Visualization
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
#        bpy.ops.object.mode_set(mode='OBJECT') # Can't assign materials in editmode

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


class ADDONNAME_OT_submit(Operator):
    bl_label = "Submit polys to potential"
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
#                        mats[i].material.diffuse_color = 0, 0, 1, 1  
                    elif enum == 'OP2':
                        face_index1.append(p.index)
#                        mat = obj.material_slots[i].material
#                        mat.diffuse_color = 1, 0, 0, 1
                    elif enum == 'OP3':
                        face_index2.append(p.index)
#                        mat = obj.material_slots[i].material
#                        mat.diffuse_color = 0, 0, 0, 1
                    elif enum == 'OP4':
                        face_index3.append(p.index)
#                        mat = obj.material_slots[i].material
#                        mat.diffuse_color = 0, 1, 0, 1
                    
        else:
            print("Object is not in edit mode.")
        surf = ADDONNAME_OT_prepare.surf
        for i in face_index:
            surf[i] = 'temp1'

        for i in face_index1:
            surf[i] = 'temp2'
        ADDONNAME_OT_submit.surf = surf

        return {'FINISHED'}

class ADDONNAME_OT_domain(Operator):
    bl_label = "Submit object to domain"
    bl_idname = "addonname.domain_operator"
    object = ''
    PONTOS_dom = []
    def execute(self, context):
        
        PONTOS_dom = []
        me = bpy.context.object.data
        bpy.ops.object.modifier_add(type='TRIANGULATE')
        bpy.ops.object.modifier_apply(modifier="Triangulate")        
        ADDONNAME_OT_domain.object = bpy.context.object
        faces = me.polygons
        verts = me.vertices
        coord = []
        elem = []
        ndom = 0
        for poly in faces:
            elem.append([poly.index+1,poly.vertices[0]+1, poly.vertices[1]+1, poly.vertices[2]+1,1])

        for v in verts:
            coord.append([v.index+1,v.co[0], v.co[1], v.co[2]])
        PONTOS_dom = coord

        ADDONNAME_OT_domain.PONTOS_dom = PONTOS_dom
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



class ADDONNAME_OT_run(Operator):
    bl_label = "Run Laplace BEM"
    bl_idname = "addonname.run_operator"

    def color_verts(my_object,T):
        vert_list = my_object.vertices
        color_map = my_object.vertex_colors.new()
        i = 0
        for poly in my_object.polygons:
            for idx in poly.loop_indices:
                loop = my_object.loops[idx]
                v = loop.vertex_index
                color_map.data[i].color = 1-T[v], 0, 0, 1
                i += 1
        return {'FINISHED'}

    def execute(self, context):
        #Solving with Liner elements
        coord = ADDONNAME_OT_prepare.coord
        elem = ADDONNAME_OT_prepare.elem
        print(f'{len(elem)} elements. Starting BEM...')
        surf = ADDONNAME_OT_submit.surf
        print(f"...{len(surf)} boundary conditions...")
        k = 1
        t_inicio = time.time()
        print('\nPrograma iniciado.')
        CCSup = {'flux0': 0., 'temp1': 1., 'temp2': 0.}
        # Cria a malha
        #malha = meshio.read(arquivo + '.msh') # Lê a malha do arquivo .msh
        # Cria a matriz NOS a partir da malha
        #NOS = malha.points
        NOS = coord
          # NOS: Matriz [NNx3] que contém as coordenadas dos vértices da malha criada,
          #      onde NN é o número de nós do problema
        # Cria a matriz ELEM a partir da malha
        #ELEM = malha.cells_dict['triangle']
        ELEM = elem
          # ELEM: Matriz [NEx3] que contém os números dos nós que formam cada elemento,
          #       onde NE é o número de elementos do problema
        print('Número de nós:',NOS.shape[0])
        print('Número de elementos:',ELEM.shape[0])
        # Cria a matriz de condições de contorno dos elementos
        CDC = contorno.gera_cdc(elem,surf,CCSup)
        #CDC = elem_bc
          # CDC: Matriz [NEx3] que contém a condição de contorno de cada nó dos
          #      elementos
        # MONTAGEM E SOLUÇÃO DO SISTEMA
        # Calcula as matrizes H e G
        npg_s = 8   # Número de pontos de Gauss para a integração singular
        npg_r = 6   # Número de pontos de Gauss para a integração regular
        qsi,w = np.polynomial.legendre.leggauss(npg_r); # Pontos e pesos de Gauss para a integração regular (triângulo)
        qsi_quad,w_quad = np.polynomial.legendre.leggauss(npg_s) # Pontos e pesos de Gauss para a integração singular (quadrilátero)
        nelem = ELEM.shape[0]           # Número de elementos
        nnos = NOS.shape[0]             # Número de nós
        H=np.zeros((nnos, nnos))        
        G=np.zeros((nnos,3*nelem))    
        t_gera_dados=time.time()-t_inicio
        t_matriz_inicio = time.time()
        print('Calculando as matrizes H e G.')
        ncores = 8      #número de threads
        H[:,:],G[:,:] = integ.cal_HeG(NOS, ELEM, k,qsi,w,qsi_quad,w_quad,ncores) # Cython
        t_matriz=time.time()-t_matriz_inicio # tempo para montagens das matrizes H e G
          # H: Matriz [NNxNN] que contém o resultado da integração de q* no contorno
          # G: Matriz [NNx3NE] que contém o resultado da integração de T* no contorno
        print('Aplicando as condições de contorno.')
        t_aplica_cdc_inicio = time.time()
        A, b, T_pr = sistema.aplica_cdc(G, H, NOS, ELEM, CDC)
        t_aplica_cdc = time.time()-t_aplica_cdc_inicio 
          # A: Matriz [NNxNN] contendo colunas de H e G
          # b: Vetor [NNx1] resultante da multiplicação de N colunas de H e G pelas
          #    CDC's conhecidas
        print('Resolvendo o sistema linear.')
        t_sistema_inicio = time.time()
        # Resolve o sistema de equações A.x = b
        x = np.linalg.solve(A, b)
          # x: Vetor [NNx1] que contém os termos calculados (antes desconhecidos) de
          #    temperatura e fluxo
        t_sistema=time.time()-t_sistema_inicio # tempo para montagens das matrizes H e G
        t_ordena_inicio=time.time()
        # Separa os valores de temperatura e fluxo
        print('Separando as variáveis.')
        T, q = sistema.monta_Teq(NOS, ELEM, CDC, x, T_pr)
          # T: Vetor [NNx1] que contém os valores de temperatura calculados
          # q: Vetor [3NEx1] que contém os valores de fluxo calculados
        t_ordena=time.time()-t_ordena_inicio
        print('Gerando o arquivo de pós-processamento.')
        t_saida_inicio=time.time()
        # Calcula os valores de temperatura no centróide do elemento
        T_centroide = contorno.calcTcentroide(T,ELEM,NOS)
        # Mostra os valores de temperatura no mapa de cor
        
        # Create material slot and material with the a color corresponding to a normalized interpolation of a color ramp.
        obj=bpy.context.object
        me = bpy.context.object.data
        ADDONNAME_OT_run.color_verts(me,T)
#        for poly in me.polygons:
#            i = poly.index
#            mat = obj.material_slots[i].material
#            mat.diffuse_color = 254*T_centroide[i]+1, 0, 0, 1

        print('Tempo para ler e gerar os dados iniciais:',t_gera_dados)
        print('Tempo para montagem das matrizes H e G:',t_matriz)
        print('Tempo para aplicas as condições de contorno:',t_aplica_cdc)
        print('Tempo para resolver o sistema linear:',t_sistema)
        print('Tempo para ordenar os dados:',t_ordena)
        #print('Tempo para gerar o arquivo de pós-processamento:',t_saida)
        print('Tempo de processamento:',time.time()-t_inicio,'s.\nPrograma finalizado.')

        return {'FINISHED'}


class ADDONNAME_OT_addcube(Operator):
    bl_label = "Add an interior cube"
    bl_idname = "addonname.addcube_operator"

    def execute(self, context):
        bpy.ops.mesh.primitive_cube_add(enter_editmode=False, align='WORLD', location=(0, 0, 0), scale=(1, 1, 1))
        bpy.ops.object.subdivision_set(level=3, relative=False)
        bpy.context.object.modifiers["Subdivision"].subdivision_type = 'SIMPLE'
        bpy.ops.object.modifier_apply(modifier="Subdivision")
#        bpy.ops.object.editmode_toggle()
#        bpy.ops.mesh.normals_make_consistent(inside=True)
#        bpy.ops.object.editmode_toggle()

        return {'FINISHED'}

class ADDONNAME_OT_addsphere(Operator):
    bl_label = "Add an exterior sphere"
    bl_idname = "addonname.addsphere_operator"

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


classes = [MyProperties, ADDONNAME_PT_main_panel, ADDONNAME_OT_submit,ADDONNAME_OT_prepare,ADDONNAME_OT_run,ADDONNAME_OT_domain,ADDONNAME_OT_addcube,ADDONNAME_OT_addsphere]
 
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
