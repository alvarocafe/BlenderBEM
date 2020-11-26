# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
using KrylovMethods, Statistics, LinearAlgebra, DelimitedFiles, PyCall, JLD
filed = "../src/waveconst3d/"
include(string(filed,"dad.jl"))
include(string(filed,"hmat.jl"))
include(string(filed,"bem_functions.jl"))
include(string(filed,"lermsh.jl"))
include(string(filed,"tree.jl"))
# Gaussian quadrature - generation of points and weights [-1,1]
npg=4; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights

files = [ "VT_A_coarsest.msh" "VT_A_coarse.msh"]
# Python - mesh.io
meshio = pyimport("meshio")
# Build the domain points
L = 140; # Length of the vocal tract
n_pint = 40; # Number of domain points
#n_pint = 10; # Number of domain points
PONTOS_int = zeros(n_pint,4);
delta = 1; # distance from both ends 
passo = (L-2*delta)/(n_pint-1);
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 delta+(i-1)*passo];
end
# Set the boundary conditions for each face. Vowel /A/ model has 30 faces
BCFace = ones(30,3);
BCFace[:,3] .= 0;
BCFace[1,:] = [1 1 1]; # Dirichlet (pressure = 1) to the Glotis
BCFace[30,:] = [30 0 0.00001]; # Dirichlet (pressure = 0) to the mouth
mshd = "./tests/data/"
t = []
for i in files
    mesh = meshio.read(string(mshd,i))
    NOS_GEO = [1:size(mesh.points,1) mesh.points]
    nelem = size(mesh.cells["triangle"],1)
    ELEM = [1:nelem mesh.cells["triangle"].+1 mesh.cell_data["triangle"]["gmsh:geometrical"]]
    CDC,NOS = gera_vars(ELEM,BCFace,NOS_GEO)
    CW = 343*1000; # Speed of sound in mm/s
    FR = 774/2/pi; # Set the frequency
    inc = [0]
    println(i," nelem = ",nelem)
    # Conventional BEM
    b1 = 1:nelem; b2 = 1:nelem;
    tmatrixc = @elapsed Ac,bc = cal_Aeb(b1,b2,[NOS,NOS_GEO,ELEM,qsi,w,CDC,FR,CW])
    tsolvec = @elapsed xc = bc\Ac
    Tc,qc = monta_Teq(CDC,xc)
    T_pintc = calc_T_pint(PONTOS_int,NOS_GEO,ELEM,Tc,qc,FR,CW,qsi,w,inc)
    println("tmatrixc = ", tmatrixc,", tsolvec = ",tsolvec)
    # H-BEM
    # max_elem = Define máximo de nós em cada folha, tal que: max_elem/2 <= nós em cada folha < max_elem
    # max_elem=floor(sqrt(2*length(NOS[:,1])))
    max_elem=10
    #println("max_elem = $max_elem")
    ttree = @elapsed Tree,child,center_row,diam,inode,ileaf = cluster(NOS[:,2:3],max_elem)
    ninterp=4 # Número de pontos de interpolação
    #η =.4 # Coeficiente relacionado a admissibilidade dos blocos
    η =.7 # Coeficiente relacionado a admissibilidade dos blocos
    allow=checa_admiss(η,center_row,diam,inode,ileaf) # mostra quais blocos são admissíveis
    block = blocks(Tree,child,allow) # Função que retorna os blocos admissiveis
    tmatrix = @elapsed hmati,bi =Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,qsi,w,CDC,FR,CW],ninterp)
    tsolve = @elapsed xi = gmres(vet->matvec(hmati,vet,block,Tree),bi,5,tol=1e-8,maxIter=1000,out=0) 
    T,q=monta_Teq(CDC,xi[1])
    T_pint = calc_T_pint(PONTOS_int,NOS_GEO,ELEM,T,q,FR,CW,qsi,w,inc)
    println("ttree = ",ttree,", tmatrix = ", tmatrix,", tsolve = ",tsolve)    
    global t = append!(t,[ttree tmatrix tsolve tmatrixc tsolvec])
    #save(string(mshd,i,".jld"),"ttree",ttree,"tmatrix",tmatrix,"tsolve","tsolve")
end # files for
