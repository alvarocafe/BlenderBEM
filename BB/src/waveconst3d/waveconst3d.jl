# Boundary element method implementation for the Helmholtz equation using constant tridimensional triangular elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Module for the parallel constant three-dimensional triangular element
# Contains the dependencies for the triangular element integration. The main function is const3D_tri.solve() which builds the influence matrices, applies the boundary conditions, solves the linear system and returns the value of the velocity potential and its flux at boundary and domain points.

module waveconst3d
using KrylovMethods, Statistics, LinearAlgebra, DelimitedFiles, PyCall, JLD, SharedArrays, DistributedArrays, Distributed
addprocs()
@everywhere filed = "./src/waveconst3d/"
@everywhere include(string(filed,"hmat.jl"))
@everywhere include(string(filed,"bem_functions.jl"))
@everywhere include(string(filed,"tree.jl"))

function solveH(info,PONTOS_int,BCFace,k,verbose=false)
    NOS_GEO,ELEM,elemint,CDC = info
    NOS = mostra_geoTRI(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    # If there's a defined BCFace, the boundary condition matrix is created. If no BCFace is defined, then the mesh reader has already built the boudnary conditions matrix.
    if isempty(BCFace) != true
	CDC = gera_CDC(ELEM,BCFace); #Monta a matriz de condicoes de contorno
    end

    # H-BEM
    # max_elem = Define máximo de nós em cada folha, tal que: max_elem/2 <= nós em cada folha < max_elem
    # max_elem=floor(sqrt(2*length(NOS[:,1])))
    max_elem=10
    #println("max_elem = $max_elem")
    ttree = @elapsed Tree,child,center_row,diam,inode,ileaf = cluster(NOS[:,2:4],max_elem)
    println(ttree)
    ninterp=4 # Número de pontos de interpolação
    #η =.4 # Coeficiente relacionado a admissibilidade dos blocos
    η =.7 # Coeficiente relacionado a admissibilidade dos blocos
    allow=checa_admiss(η,center_row,diam,inode,ileaf) # mostra quais blocos são admissíveis
    block = blocks(Tree,child,allow) # Função que retorna os blocos admissiveis
    Parallel_tmatrix = @elapsed Parallel_hmati,parallel_bi = parallel_Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,qsi,w,CDC,FR,CW],ninterp)
    println("Paralelo = ",Parallel_tmatrix)
    tmatrix = @elapsed hmati,bi =Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,qsi,w,CDC,FR,CW],ninterp)
    println("não paralelo = ",tmatrix)
    tsolve = @elapsed xi = gmres(vet->matvec(hmati,vet,block,Tree),bi,5,tol=1e-8,maxIter=1000,out=0)
    println(tsolve)
    T,q=monta_Teq(CDC,xi[1])
    T_pint = calc_T_pint(PONTOS_int,NOS_GEO,ELEM,T,q,FR,CW,qsi,w,inc)
    println("ttree = ",ttree,", tmatrix = ", tmatrix,", tsolve = ",tsolve)    

return T,q,T_pint
end

end
