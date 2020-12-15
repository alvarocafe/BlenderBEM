# -*- coding: utf-8 -*-
"""
Programa:  sistema.py
Descrição: Contém as funções para a solução do problema.
Funções:   cal_HeG
           calcula_Gs
           calcula_HeGns
           calc_vetnormal
           calc_jacobiano
           calc_fforma_quad
           calc_dfforma_quad
           calc_jacobiano_quad

@author: Gustavo Gontijo (ggontijo@gmail.com)
"""
import numpy as np
from numpy import linalg as LA

#%%  CÁLCULO DOS ELEMENTOS DAS MATRIZES H E G
def cal_HeG(NOS,ELEM,k):
    nelem = ELEM.shape[0]           # Número de elementos
    nnos = NOS.shape[0]             # Número de nós
    
    G = np.zeros((nnos, 3*nelem))   # Inicialização da matriz G
    H = np.zeros((nnos, nnos))      # Inicialização da matriz H
        
    #%% = = = = = PONTOS E PESOS DE GAUSS = = = = =
    
    npg_s = 8   # Número de pontos de Gauss para a integração singular
    npg_r = 6   # Número de pontos de Gauss para a integração regular

    qsil,w1 = np.polynomial.legendre.leggauss(npg_r); # Pontos e pesos de Gauss para a integração regular (triângulo)
    qsil = 0.5*(qsil + 1)                             # Converte os pontos de Gauss para o intervalo [0, 1]
    w1 = w1*0.5                                       # Converte os pesos de Gauss para o intervalo [0, 1]
    eta = qsil                                        # Utiliza os mesmos pontos de Gauss para as duas coordenadas
    w2 = w1                                           # Utiliza os mesmos pesos de Gauss para as duas coordenadas

    qsi_quad,w_quad = np.polynomial.legendre.leggauss(npg_s) # Pontos e pesos de Gauss para a integração singular (quadrilátero)
    
    #%% = = = = = FUNÇÕES DE FORMA = = = = =
    
    N1 = np.zeros((npg_r,npg_r))
    N2 = np.zeros((npg_r,npg_r))
    for l in range(0,npg_r):
        for m in range(0,npg_r):
            N1[l,m] = (1-eta[m])*eta[l] # qsi escrito como qsi_linha e eta
            N2[l,m] = eta[m]
    N3 = 1-N1-N2
    dNdqsi = np.array([1,0,-1])
    dNdeta = np.array([0,1,-1])
        
    #%% = = = = = CÁLCULO DOS ELEMENTOS DAS MATRIZES H e G = = = = =
    
    for pc in range(0,nelem):   # Laço sobre os elementos (pontos campo)
        nos = ELEM[pc]      # Nós que compõem o elemento
        X1 = NOS[nos[0]]    # Coordenadas do nó 1 do elemento
        X2 = NOS[nos[1]]    # Coordenadas do nó 2 do elemento
        X3 = NOS[nos[2]]    # Coordenadas do nó 3 do elemento
    
        n = calc_vetnormal(X1,X2,X3)    # Vetor unitário normal ao elemento
        J = calc_jacobiano(X1,X2,X3,dNdqsi,dNdeta)  # Jacobiano do elemento
        
        for pf in range(0,nnos): # Laço sobre os pontos fonte
            Xd = NOS[pf]         # Coordenadas do ponto fonte
            
            if pf in nos: # Integração singular (o ponto fonte pertence ao elemento)
                g = calcula_Gs(X1,X2,X3,Xd,qsi_quad,w_quad,k)
                h = np.zeros((g.shape[0]))
            else: # Integração regular (o ponto fonte não pertence ao elemento)
                g,h = calcula_HeGns(X1,X2,X3,Xd,eta,w2,n,k,J,N1,N2,N3)
            
            for nolocal in range(0,3):
                noglobal = ELEM[pc,nolocal]   # Índice da matriz global H
                H[pf,noglobal] = H[pf,noglobal] + h[nolocal]

                
            G[pf,3*pc:3*pc+3] = g

            
    #%% = = = = = DIAGONAL DA MATRIZ H = = = = =
    # Calcula os termos da diagonal da matriz H (consideração de corpo a temperatura constante)
    for m in range(0,nnos):
       H[m,m] = 0             # Zera a diagonal principal
       for n in range(0,nnos):
          if n != m:
             H[m,m] = H[m,m] - H[m,n]
             
    print(H[2,:])             
    
    return H, G


#%%  CÁLCULO DO VETOR NORMAL
def calc_vetnormal(X1,X2,X3):
    # Function que calcula o vetor unitário normal ao elemento

    v1 = X3 - X2           # Vetor formado pela aresta 3-2 do elemento
    v2 = X1 - X2           # Vetor formado pela aresta 1-2 do elemento
    n = np.cross(v1, v2)   # Produto vetorial entre v1 e v2 (vetor normal ao elemento)
    n = n/LA.norm(n)         # Vetor unitário normal ao elemento

    return n


#%%  CÁLCULO DO JACOBIANO
def calc_jacobiano(X1,X2,X3,dNdqsi,dNdeta):
    # Cálculo do Jacobiano para elementos triangulares lineares
    # Baseado na página 221 do livro do Dominguez
    
    dxdqsi = X1[0]*dNdqsi[0] + X2[0]*dNdqsi[1] + X3[0]*dNdqsi[2]
    dydqsi = X1[1]*dNdqsi[0] + X2[1]*dNdqsi[1] + X3[1]*dNdqsi[2]
    dzdqsi = X1[2]*dNdqsi[0] + X2[2]*dNdqsi[1] + X3[2]*dNdqsi[2]

    dxdeta = X1[0]*dNdeta[0] + X2[0]*dNdeta[1] + X3[0]*dNdeta[2]
    dydeta = X1[1]*dNdeta[0] + X2[1]*dNdeta[1] + X3[1]*dNdeta[2]
    dzdeta = X1[2]*dNdeta[0] + X2[2]*dNdeta[1] + X3[2]*dNdeta[2]

    g1 = dydqsi*dzdeta - dzdqsi*dydeta
    g2 = dzdqsi*dxdeta - dxdqsi*dzdeta
    g3 = dxdqsi*dydeta - dydqsi*dxdeta
    J = np.sqrt(g1**2.0 + g2**2.0 + g3**2.0)
   
    return J


#%%  INTEGRAÇÃO SINGULAR DA MATRIZ G
def calcula_Gs(X1,X2,X3,Xd,qsi_quad,w_quad,k):
    """
    Function:  calcula_Gs
    Descrição: Integração singular da matriz G. O elemento triangular é
               transformado em um elemento quadrilateral degenerado. Os dois 
               primeiros nós deste quadrilátero são coincidentes (formam um 
               só vértice do triângulo) e correspondem ao ponto onde existe a
               singularidade, ou seja, ao ponto fonte. Isto faz com que haja 
               uma concentração de pontos de integração junto à 
               singularidade, além do jacobiano ser igual a zero na 
               singularidade. No caso da matriz H, os elementos da diagonal
               são calculados pela consideração de corpo a temperatura 
               constante.
    Autor:     Gustavo Gontijo, adaptado de Éder Lima de Albuquerque
       
    Última modificação: 05/05/2014 - 19h23min
    """   
    npg = qsi_quad.shape[0]   # Número de pontos de integração
    g = np.zeros(3)           # Inicialização da matriz g
    J = np.zeros(npg)
    G=np.zeros((npg,npg))
    
    # Transformação do triângulo em quadrilátero
    
    X1t = Xd   # coordenadas do primeiro nó do quadrilatero degenerado
    X4t = Xd   # coordenadas do quarto nó do quadrilatero degenerado
    
    if np.all(X1 == Xd):   # O ponto fonte está no nó 1
        X2t = X2             # coordenadas do segundo nó do quadrilatero degenerado
        X3t = X3             # coordenadas do terceiro nó do quadrilatero degenerado
        pf = 1               # Sinaliza que o ponto fonte está sobre o primeiro nó
    elif np.all(X2 == Xd): # O ponto fonte está no nó 2
        X2t = X3             # coordenadas do segundo nó do quadrilatero degenerado
        X3t = X1             # coordenadas do terceiro nó do quadrilatero degenerado
        pf = 2               # Sinaliza que o ponto fonte está sobre o segundo nó
    elif np.all(X3 == Xd): # O ponto fonte está no nó 3
        X2t = X1             # coordenada do segundo nó do quadrilatero degenerado
        X3t = X2             # coordenada do terceiro nó do quadrilatero degenerado
        pf = 3               # Sinaliza que o ponto fonte está sobre o terceiro nó

    # Integração regular do elemento
    for l in range(0,npg):       # Laço sobre a primeira variável de integração
        for m in range(0,npg):   # Laço sobre a segunda variável de integração
            # Cálculo das funções de forma
            N = calc_fforma_quad(qsi_quad[l],qsi_quad[m])
            # Coordenadas do ponto campo
            Xc = N[0]*X1t + N[1]*X2t + N[2]*X3t + N[3]*X4t

            
            # Jacobiano (varia ao longo do elemento degenerado)
            if m == 0:
                J[l] = calc_jacobiano_quad(X1t,X2t,X3t,X4t,qsi_quad[l],qsi_quad[m])

            # Solução fundamental: início
            R = Xc - Xd
            r = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)
            Tast = 1.0/(4.0*k*np.pi*r)
            # Solução fundamental: fim

            # Integral da matriz G
            if pf==1:
                g = g + Tast*np.array([N[0]+N[3],N[1],N[2]])*w_quad[l]*w_quad[m]*J[l]
                G[l,m]=(N[0]+N[3])*w_quad[l]*w_quad[m]*J[l]*Tast
            elif pf==2:
                g = g + Tast*np.array([N[2],N[0]+N[3],N[1]])*w_quad[l]*w_quad[m]*J[l]
            else:
                g = g + Tast*np.array([N[1],N[2],N[0]+N[3]])*w_quad[l]*w_quad[m]*J[l]
    return g


#%%  FUNÇÕES DE FORMA DO ELEMENTO QUADRILATERAL
def calc_fforma_quad(qsi,eta):
    """
    Dados de entrada:
    qsi, eta - pontos de Gauss onde as funções de forma são calculadas.
    Dados de saída:
    [N] - Funções de forma para um elemento quadrilateral linear calculadas em (qsi, eta).
    """
    N = (1./4.)*np.array([(1. - qsi)*(1. - eta),
                          (1. + qsi)*(1. - eta),
                          (1. + qsi)*(1. + eta),
                          (1. - qsi)*(1. + eta)])
    
    return N


#%%  JACOBIANO DO ELEMENTO QUADRILATERAL
def calc_jacobiano_quad(X1,X2,X3,X4,qsi,eta):
   dNdqsi, dNdeta = calc_dfforma_quad(qsi,eta) # Calcula a derivada das funções de forma
   
   dxdqsi = X1[0]*dNdqsi[0] + X2[0]*dNdqsi[1] + X3[0]*dNdqsi[2] + X4[0]*dNdqsi[3]
   dydqsi = X1[1]*dNdqsi[0] + X2[1]*dNdqsi[1] + X3[1]*dNdqsi[2] + X4[1]*dNdqsi[3]
   dzdqsi = X1[2]*dNdqsi[0] + X2[2]*dNdqsi[1] + X3[2]*dNdqsi[2] + X4[2]*dNdqsi[3]
   
   dxdeta = X1[0]*dNdeta[0] + X2[0]*dNdeta[1] + X3[0]*dNdeta[2] + X4[0]*dNdeta[3]
   dydeta = X1[1]*dNdeta[0] + X2[1]*dNdeta[1] + X3[1]*dNdeta[2] + X4[1]*dNdeta[3]
   dzdeta = X1[2]*dNdeta[0] + X2[2]*dNdeta[1] + X3[2]*dNdeta[2] + X4[2]*dNdeta[3]

   g1 = dydqsi*dzdeta - dzdqsi*dydeta
   g2 = dzdqsi*dxdeta - dxdqsi*dzdeta
   g3 = dxdqsi*dydeta - dydqsi*dxdeta
   J = np.sqrt(g1**2.0 + g2**2.0 + g3**2.0)

   return J


#%%  DERIVADAS DAS FUNÇÕES DE FORMA DO ELEMENTO QUADRILATERAL
def calc_dfforma_quad(qsi,eta):
   dNdqsi = (1./4.)*np.array([-(1. - eta),
                               (1. - eta),
                               (1. + eta),
                              -(1. + eta)])

   dNdeta = (1./4.)*np.array([-(1. - qsi),
                              -(1. + qsi),
                               (1. + qsi),
                               (1. - qsi)])
   
   return dNdqsi, dNdeta


#%%  INTEGRAÇÃO REGULAR DAS MATRIZES H e G
#def calcula_HeGns(X1,X2,X3,Xd,qsil,w1,eta,w2,n,k,J,N1,N2,N3):
def calcula_HeGns(X1,X2,X3,Xd,eta,w2,n,k,J,N1,N2,N3):
   """
   Integração numérica não-singular de elementos de contorno triangulares
   lineares contínuos
   """
   g = np.array([0, 0, 0])    # Inicializa o somatório de g
   h = np.array([0, 0, 0])    # Inicializa o somatório de h

#   n_pint1 = qsil.shape[0]    # Nro de pontos de integração na direção qsi
   n_pint = eta.shape[0]     # Nro de pontos de integração na direção eta

   for l in range(0,n_pint):     # Laço sobre os pontos de integração
      for m in range(0,n_pint):  # Laço sobre os pontos de integração
         Xc = N1[l,m]*X1 + N2[l,m]*X2 + N3[l,m]*X3 # coordenadas dos pontos de integração
        
         # Solução fundamental: início
         R = Xc - Xd
         r = np.sqrt(R[0]**2 + R[1]**2 + R[2]**2)
         Tast = 1.0/(4.0*k*np.pi*r)
         qast = (R[0]*n[0] + R[1]*n[1] + R[2]*n[2])/(4.0*np.pi*r**3.0)
         # Solução fundamental: fim
        
         # Integral da matriz H
         h = h + qast*np.array([N1[l,m],N2[l,m],N3[l,m]])*(1-eta[m])*w2[l]*w2[m]*J
         # Integral da matriz G
         g = g + Tast*np.array([N1[l,m],N2[l,m],N3[l,m]])*(1-eta[m])*w2[l]*w2[m]*J

   return g, h


#%%  APLICAÇÃO DAS CONDIÇÕES DE CONTORNO
def aplica_cdc(G,H,NOS,ELEM,CDC):
   """
   Programa:  aplica_cdc.m
   Descrição: Aplica as condições de contorno e separa os valores conhecidos
              dos desconhecidos, gerando a matriz [A] (que contém os
              coeficientes do vetor de valores desconhecidos {x}) e o vetor
              {b} (que contém o resultado da multiplicação da matriz dos
              coeficientes dos valores conhecidos pelo vetor de valores
              conhecidos). Cria também o vetor T_pr, que tem marcado com o
              valor 1 as posições do vetor cujo índice é o número de um nó
              que tem temperatura prescrita.
   Autor:     Gustavo Gontijo

   Última modificação: 07/03/2014 - 20h19min
   """
   nnos = NOS.shape[0]    # Nro de nós
   nelem = ELEM.shape[0]  # Nro de elementos
   todosvalores = np.zeros(3*nelem)
   T_pr = np.zeros((nnos,nelem), dtype=bool)  # Matriz contendo valores "False"
   
   for i in range(0,nelem):  # For sobre os elemento
       #Lê a condição de contorno deste elemento
      tipo_cdc = CDC[i,0]
      valor_cdc = CDC[i,1]
    
      indH = ELEM[i] # Índices das colunas da matriz H que serão trocadas

      indG = (3*i) + np.array([0, 1, 2]) # Índice das colunas da matriz G que serão trocadas
    
      if tipo_cdc == 0:    # Se o valor prescrito é de potencial
         # Então a coluna da matriz H de todos os nós deste elemento deve ser trocada
        
         # Atualiza o vetor T_pr
         nos = ELEM[i]
         T_pr[nos[0],i] = 1
         T_pr[nos[1],i] = 1
         T_pr[nos[2],i] = 1
         # Se o nó 'x' tem temperatura prescrita pelo elemento 'el', então,
         # a matriz T_pr, na posição (x,el), recebe o valor 1. Ao final, se 
         # o nó x não tem temperatura prescrita em nenhum elemento, o 
         # conteúdo de T_pr(x,:) é zero.
        
         # Verifica se o nó já teve o potencial prescrito por um elemento
         # anterior. Caso afirmativo, não realiza a troca de colunas, mas
         # sua soma.
        
         ja_trocado = np.zeros(3)

         if i != 0:
             for j_no in range(0,3):
                 for j_el in range(0,i):   # For sobre os elementos anteriores ao atual
                     if T_pr[nos[j_no],j_el] == 1:
                         ja_trocado[j_no] = 1
        
         # Procede a troca das colunas
         troca=np.zeros(nnos)
         for j in range(0,3):
             if ja_trocado[j] == 0:   # Esta é a primeira vez que o potencial
                                      # é prescrito para este nó
                 troca[:] = G[:,indG[j]]
                 G[:,indG[j]] = -H[:,indH[j]]
                 H[:,indH[j]] = -troca
             else:                  # O potencial já foi prescrito para este nó
                 H[:,indH[j]] = H[:,indH[j]] - G[:,indG[j]]
                 G[:,indG[j]] = 0

      # Coloca o valor da condição de contorno deste elemento no vetor
      # todosvalores, na posição dos nós correspondentes.
      for j in range(0,3):
          todosvalores[indG[j]] = valor_cdc

   b = np.dot(G, todosvalores)
   
   return H, b, T_pr
   

#%%  SEPARAÇÃO DOS VALORES DE T e Q
def monta_Teq(NOS, ELEM, CDC, x, T_pr):
   """
   Programa:  monta_Teq.m
   Descrição: Separa o vetor de valores calculados {x} e os valores
              prescritos [CDC] em valores de potencial {T} e valores de
              fluxo {q}.
   Autor:     Gustavo Gontijo

   Última modificação: 12/03/2014 - 18h10min
   """
   
   nnos = NOS.shape[0]     # Nro de nós
   nelem = ELEM.shape[0]   # Nro de elementos
   T = np.zeros(nnos)      # Inicia o vetor T (potenciais)
   q = np.zeros(3*nelem)   # Inicia o vetor q (fluxos)
   
   # Coloca os valores prescritos nos vetores T e q
   for i in range(0,nelem):
       tipo_cdc = CDC[i,0]
       valor_cdc = CDC[i,1]
       if tipo_cdc == 0:    # O potencial é prescrito
           indT = ELEM[i]
           for k in range(0,3):
               T[indT[k]] = valor_cdc
       else:                # O fluxo é prescrito
           indq = (i*3) + np.array([0, 1, 2])
           for k in range(0,3):
               q[indq[k]] = valor_cdc

   # Adiciona nos vetores T e q os valores calculados (estão em {x})
   for i in range(0,nnos):
      sw = 0
      for j in range(0,nelem):
         if T_pr[i,j] == 1:  # O nó teve potencial prescrito no elemento j,
                             # portanto, o valor calculado é de fluxo.
            sw = 1
            for k in range(0,3):
               if ELEM[j,k] == i:
                  indq = (j*3) + k
            
            q[indq] = x[i]
   
      if sw == 0:  # O nó não teve potencial prescrito em nenhum elemento,
                   # portanto, o valor calculado é de potencial.
        T[i] = x[i]
   
   return T, q