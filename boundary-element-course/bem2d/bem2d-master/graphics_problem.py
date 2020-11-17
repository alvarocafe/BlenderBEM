# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 14:29:44 2017

@author: Éder
"""
import matplotlib.tri as triang
import matplotlib.pyplot as plt
from matplotlib import cm
import geometry
import numpy as np

# def show_problem(node_bound,node_int,node_med,normal,elem,coord,bcs,T,q,Tint,tri):
def show_problem(node_med,normal,coord,bcs,tri):
    plt.close("all")
    nelem=bcs.shape[0]
#    node_ends=coord[node_bound]
    indknownT=[]
    indknownq=[]
    for i in range(nelem):
        if(bcs[i,0]==0):
            indknownT.append(i)
        else:
            indknownq.append(i)
    xmax=np.max(coord[:,0])
    ymax=np.max(coord[:,1])
    xmin=np.min(coord[:,0])
    ymin=np.min(coord[:,1])
    deltax=xmax-xmin
    deltay=ymax-ymin
    ax =plt.gca() # get current axes
    plt.grid("on")
    plt.quiver(node_med[indknownq,0],node_med[indknownq,1],normal[indknownq,0]*bcs[indknownq,1],
               normal[indknownq,1]*bcs[indknownq,1],color="blue",width=0.002,
                     scale=100, headaxislength=0)
    plt.quiver(node_med[indknownT,0],node_med[indknownT,1],normal[indknownT,0]*bcs[indknownT,1],
               normal[indknownT,1]*bcs[indknownT,1],color="red",width=0.002,
                     scale=100, headaxislength=0)
    plt.plot(node_med[indknownT,0],node_med[indknownT,1],"ro",markersize=4)	# Plot the node of the elements
    plt.plot(node_med[indknownq,0],node_med[indknownq,1],"bo",markersize=4)	# Plot the node of the elements
    plt.triplot(coord[:,0], coord[:,1], tri, color=(0.0,0.,0.),linewidth=0.4)
    ax.set_xlim(xmin-.15*deltax,xmax+.15*deltax)
    ax.set_ylim(ymin-.15*deltay,ymax+.15*deltay)
    
    
def show_results(node_all,node_int,node_med,elem,coord,T,q,Tint):
    elem_local=compute_bounds(elem,node_all)
    coord_tri2=np.concatenate((node_med[:,0:2], coord[node_int,0:2]), axis=0)
    trimesh = triang.Triangulation(coord_tri2[:,0],coord_tri2[:,1])
    tri2=trimesh.triangles    
    plt.figure()
    cor=np.concatenate((T,Tint))
    centroid=(coord_tri2[tri2[:,0],:]+coord_tri2[tri2[:,1],:]+coord_tri2[tri2[:,2],:])/3.0
    ind = geometry.inpoly(centroid[:,0:2],node_med[:,0:2], elem_local)
    tri3 = tri2[ind,:]
    plt.triplot(coord_tri2[:,0], coord_tri2[:,1], tri3, color=(0.0,0.,0.),linewidth=0.4)
    plt.tricontourf(coord_tri2[:,0], coord_tri2[:,1], tri3, cor,cmap=cm.jet)
    plt.axis("equal")
    plt.title('Temperature')
    plt.colorbar()
    plt.show()
    
def compute_bounds(elem,nodes_all):
    nelem=elem.shape[0]
    elem_local=np.zeros((nelem,2),dtype=np.int)
    nnodes=len(nodes_all)
    nodes_local=np.zeros(nnodes)
    seq_nodes=elem[:,0]
    nnodes_local=len(seq_nodes)
    for t in range(0,nnodes_local):
        inode_global=seq_nodes[t]
        nodes_local[inode_global]=t
    for t in range(0,nelem):
        inode1=elem[t,0]
        inode2=elem[t,1]
        inode1_local=nodes_local[inode1]
        inode2_local=nodes_local[inode2]
        elem_local[t]=[inode1_local,inode2_local]
    return elem_local    

