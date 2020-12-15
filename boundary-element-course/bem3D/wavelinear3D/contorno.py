import numpy as np

def gera_CDC(CCSup, SUP, ELEM):
    CDC = np.zeros((len(ELEM), 4))
    
    # Estabelece fluxo igual a zero para todos os elementos
    CDC[:,0] = 1
    
    for superf, tipo_valor in CCSup.items():
       el = np.where(SUP == superf)
       CDC[el,0] = tipo_valor[0]
       CDC[el,1:4] = tipo_valor[1]

    return CDC

def gera_elem_cdc(malha, CCSup):
    superf = malha.cell_data_dict["gmsh:physical"]["triangle"]
    lista_cdc = malha.field_data
    tipo_cdc = []
    valor_cdc = []
    bc_seg = []
    for i in lista_cdc:
        valor_cdc.append(CCSup[i])
        bc_seg.append(malha.field_data[i][0])
        if(i[0:4] == 'temp'):
            tipo_cdc.append(0)
        else:
            tipo_cdc.append(1)
    nn = len(superf)
    elem_bc = np.zeros([nn,4])
    for i in range(0,nn):
        isurf = superf[i]
        elem_bc[i,0] = tipo_cdc[bc_seg.index(isurf)]
        elem_bc[i,1:4] = valor_cdc[bc_seg.index(isurf)]
    return elem_bc

def gera_cdc(superf,lista_cdc, CCSup):
#    superf = malha.cell_data_dict["gmsh:physical"]["triangle"]
#    lista_cc = malha.field_data
    tipo_cdc = []
    valor_cdc = []
    bc_seg = []
    for i in range(len(lista_cdc)):
        valor_cdc.append(CCSup[lista_cdc[i]])
        bc_seg.append(superf[i,:])
        if(lista_cdc[i][0:4] == 'temp'):
            tipo_cdc.append(0)
        else:
            tipo_cdc.append(1)
    nn = len(superf)
    elem_bc = np.zeros([nn,4])
    for i in range(0,nn):
        isurf = superf[i]
        elem_bc[i,0] = tipo_cdc[i]
        elem_bc[i,1:4] = valor_cdc[i]
    return elem_bc
    
def calcTcentroide(T, ELEM, NOS):
   """
   Calcula o potencial no centroide do elemento através de uma média ponderada 
   pelo inverso da distância do nó até o centróide.
   """
   n_elem = ELEM.shape[0]
   T_centroide = np.zeros(n_elem)
   for i in range(0,n_elem):
      # Número dos nós que compõem o elemento
      no1 = ELEM[i,0]
      no2 = ELEM[i,1]
      no3 = ELEM[i,2]
      
      # Calcula o centróide do elemento
      xc = (NOS[no1,0] + NOS[no2,0] + NOS[no3,0])/3
      yc = (NOS[no1,1] + NOS[no2,1] + NOS[no3,1])/3
      zc = (NOS[no1,2] + NOS[no2,2] + NOS[no3,2])/3
      
      r1 = np.sqrt((NOS[no1,0]-xc)**2 + (NOS[no1,1]-yc)**2 + (NOS[no1,2]-zc)**2)
      r2 = np.sqrt((NOS[no2,0]-xc)**2 + (NOS[no2,1]-yc)**2 + (NOS[no2,2]-zc)**2)
      r3 = np.sqrt((NOS[no3,0]-xc)**2 + (NOS[no3,1]-yc)**2 + (NOS[no3,2]-zc)**2)
      
      T_centroide[i] = (T[no1]/r1 + T[no2]/r2 + T[no3]/r3)/(1/r1 + 1/r2 + 1/r3)
   
   return T_centroide