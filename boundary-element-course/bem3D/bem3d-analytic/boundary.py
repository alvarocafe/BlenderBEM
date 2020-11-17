"""

"""
import numpy as np

def elem(surf_bc, surf, elem):
    elem_bc = np.zeros((len(elem), 2))
    #all elements with fluxo 0
    elem_bc[:,0] = 1
    for bc_surf,  bc_value  in surf_bc.items():
        ix = np.where(surf == bc_surf)[0]
        elem_bc[ix] = bc_value 

    return elem_bc

def correct_T(T,extremes):
    for i in range(T.shape[0]):
        if(T[i]>extremes[1]):
            T[i]=extremes[1]
        elif(T[i]<extremes[0]):
            T[i]=extremes[0]
    return T

    
    
# usinfg a for loop
#        for i, elem_surf in enumerate(surf):
#            if elem_surf == bc_surf:
#                elem_bc[i] = [bc_type, bc_value]