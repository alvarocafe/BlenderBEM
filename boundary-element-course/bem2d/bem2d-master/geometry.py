# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 09:27:11 2016

@author: mansour
"""
import numpy as np
from mesh import Mesh


def comp_node_and_normal(elem, nodes):
    # element half point - list
    nelem=elem.shape[0]
    node_med = np.zeros((nelem,2))
    normal = np.zeros((nelem,2))
    # loop through elements
    for i in range(0,nelem):
        node1=elem[i,0]
        node2=elem[i,1]
            # first node
        x1=nodes[node1,0]
        y1=nodes[node1,1]
        #    print(x1,y1)
            # second node
        x2=nodes[node2,0]
        y2=nodes[node2,1]
        #    print(x2,y2)
            # font nodes
        node_med[i,]=([(x2+x1)/2, (y2+y1)/2])
        L=np.sqrt((x2-x1)**2+(y2-y1)**2)
        sx=(x2-x1)/L
        sy=(y2-y1)/L
        normal[i,:]=[sy,-sx]        
            # node_med.append([(x2+x1)/2, (y2+y1)/2])    
            # print('element half points', node_med[line][0])    
    return node_med,normal

def compute_inodes(file):
    mymesh = Mesh() #criou objeto: instanciou a classe
    #chamei o metodo do objeto
    mymesh.read_msh(file + '.msh')
    coord = mymesh.Verts
    elem = mymesh.Elmts[1][1]-1
    segmentos = mymesh.Elmts[1][0]-1
    tri = mymesh.Elmts[2][1]-1
    inode_bound = np.unique(elem) # index of boundary nodes
    inode_all = np.unique(tri) # index of all nodes (boundary and interiors)
    inode_int = [inode for inode in inode_all if inode not in inode_bound]
    return inode_bound,inode_int,inode_all,coord,elem,segmentos,tri

def inpoly(bb,node,edge):
    reltol = 1.0e-12

# #  PRE-PROCESSING
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    nn  = bb.shape[0]
    nc = edge.shape[0]

#  Choose the direction with the biggest range as the "y-coordinate" for the
#  test. This should ensure that the sorting is done along the best
#  direction for long and skinny problems wrt either the x or y axes.
    dxy = np.max(bb,0)-np.min(bb,0);
    if dxy[0]>dxy[1]:
   #  Flip co-ords if x range is bigger
        bb = bb[:,[1,0]]
        node = node[:,[1,0]]


#  Polygon bounding-box
    dxy = np.max(node,0)-np.min(node,0)
    tol = reltol*np.min(dxy)
    if tol==0.0:
        tol = reltol

#  Sort test points by y-value
    if(bb.shape[0]>1):
        i = np.argsort(bb[:,1])
        y= bb[i,1]
        x = bb[i,0]
    else:
        y=np.array([bb[0,1]])
        x=np.array([bb[0,0]])
        i=0
# #  MAIN LOOP
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    vet = np.zeros(nn, dtype=bool) 
    cn = np.zeros(nn, dtype=bool)      #  Because we're dealing with mod(cn,2) we don't have
                  #  to actually increment the crossing number, we can
                     #  just flip a logical at each intersection (faster!)
    on = cn[:]
    for k in range(0,nc):         #  Loop through edges
    #  Nodes in current edge
        n1 = edge[k,0]
        n2 = edge[k,1]

    #  Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
    #            - also get xmin = min(x1,x2), xmax = max(x1,x2)
        y1 = node[n1,1]
        y2 = node[n2,1]
        if y1<y2:
            x1 = node[n1,0]
            x2 = node[n2,0]
        else:
            yt = y1
            y1 = y2
            y2 = yt
            x1 = node[n2,0]
            x2 = node[n1,0]
        if x1>x2:
            xmin = x2
            xmax = x1
        else:
            xmin = x1
            xmax = x2
    #  Binary search to find first point with y<=y1 for current edge
        if y[0]>=y1:
            start = 0
        elif y[nn-1]<y1:
            start = nn  
        else:
            lower = 0
            upper = nn-1
            for j in range(0,nn):
                start = round(1/2*(lower+upper))
                if y[start]<y1:
                    lower = start+1
                elif y[start-1]<y1:
                    break
                else:
                    upper = start-1
    #  Loop through points
        for jj in range(start,nn):
        #  Check the bounding-box for the edge before doing the intersection
        #  test. Take shortcuts wherever possible!
            Y = y[jj]   #  Do the array look-up once & make a temp scalar
            if Y<=y2:
                X = x[jj]   #  Do the array look-up once & make a temp scalar
                if X>=xmin:
                    if X<=xmax:

               #  Check if we're "on" the edge
                        on[jj] = on[jj] or (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<=tol)

               #  Do the actual intersection test
                        if ((Y<y2) and (y2-y1)*(X-x1)<(Y-y1)*(x2-x1)):
                            cn[jj] = ~cn[jj]
                elif Y<y2:   #  Deal with points exactly at vertices
               #  Has to cross edge
                   cn[jj] = ~cn[jj]
            else:
               #  Due to the sorting, no points with >y
               #  value need to be checked
               break
#  Re-index to undo the sorting
    vet[i] = np.logical_or(cn, on)

    return vet