# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
#module const3D
using KrylovMethods, Statistics, LinearAlgebra, DelimitedFiles, PyCall, JLD
include("dad.jl")
include("hmat.jl")
include("bem_functions.jl")
include("lermsh.jl")
include("tree.jl")
# Gaussian quadrature - generation of points and weights [-1,1]
npg=4; # Number of integration points
qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights

function solve(info,PONTOS_int,fc,BCFace,k)
    ## CBIE - Conventional Boundary Integral Equation
    NOS_GEO,NOS,ELEM,CDC = info
    nnos = size(NOS,1)  # Number of physical nodes, same as elements
    b1 = 1:nnos # Array containing all the indexes for nodes and elements
    npg=6
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    x = A\b # Solves the linear system
    phi,qphi = monta_phieq(CDC,x) # Applies the boundary conditions
    phi_dom = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,fc,qsi,w,k)
    return phi, qphi, phi_dom, phi_dom
end
function solveH(info,PONTOS_int,fc,BCFace,k)
    ## H-Matrix BEM - Interpolation using Lagrange polynomial
    NOS_GEO,NOS,ELEM,CDC = info
    Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
    npg=6
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    vet->matvec(A,vet,block,Tree)
    x = gmres(vet,b)
    phi,qphi = monta_phieq(CDC,x[1]) # Applies the boundary conditions
    phi_dom = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,fc,qsi,w,k)
    return phi,qphi,phi_dom,phi_dom
end

files = ["VT_A_coarse.msh" "VT_A_coarsest.msh" "VT_A_fine.msh" "VT_A_very_fine.msh" "VT_A.msh"]
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
BCFace[1,:] = [1 0 1]; # Dirichlet (pressure = 1) to the Glotis
BCFace[30,:] = [30 0 0]; # Dirichlet (pressure = 0) to the mouth
mshd = "tests/data/"
t = []
for i in files
    mesh = meshio.read(string(mshd,i))
    NOS_GEO = [1:size(mesh.points,1) mesh.points]
    ELEM = [1:size(mesh.cells["triangle"],1) mesh.cells["triangle"].+1 mesh.cell_data["triangle"]["gmsh:geometrical"]]
    CDC,NOS = gera_vars(ELEM,BCFace,NOS_GEO)
    FR = 1;CW = 1;inc = [0]
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
    save(string(mshd,i,".jld"),"ttree",ttree,"tmatrix",tmatrix,"tsolve","tsolve")
end # files for
#end # module const3D
