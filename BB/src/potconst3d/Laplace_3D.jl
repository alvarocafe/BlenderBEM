using KrylovMethods
include("dad.jl")
include("hmat.jl")
include("bem_functions.jl")
#include("pre_proc.jl")
include("tree.jl")
using Statistics
using LinearAlgebra


NOS_GEO,ELEM,CCFace,k=dad3() # Arquivo de entrada de dados
# dadGID # Arquivo de entrada de dados

nelem=size(ELEM,1)
println("Nũmero de nós: $nelem")


CDC,NOS = gera_vars(ELEM,CCFace,NOS_GEO);   # Gera a matriz de condições de contorno
# Mostra a geometria do problema

npg=4;
qsi,w = Gauss_Legendre(-1,1,npg) # Gera os pontos e pesos de Gauss

# max_elem = Define máximo de nós em cada folha, tal que: max_elem/2 <= nós em cada folha < max_elem
# max_elem=floor(sqrt(2*length(NOS[:,1])))
max_elem=10
println("max_elem = $max_elem")
Tree,child,center_row,diam,inode,ileaf = cluster(NOS[:,2:3],max_elem)

ninterp=4 # Número de pontos de interpolação
#η =.4 # Coeficiente relacionado a admissibilidade dos blocos
η =.7 # Coeficiente relacionado a admissibilidade dos blocos

allow=checa_admiss(η,center_row,diam,inode,ileaf) # mostra quais blocos são admissíveis

block = blocks(Tree,child,allow) # Função que retorna os blocos admissiveis


hmati,bi =Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,qsi,w,CDC,k],ninterp)

xi = gmres(vet->matvec(hmati,vet,block,Tree),bi,5,tol=1e-8,maxIter=1000,out=0) 

T,q=monta_Teq(CDC,xi[1])

A=montacheia(hmati,block,Tree,nelem)
# L,HH,G2=cal_Aeb_interp(Tree[14],Tree[8],[NOS,NOS_GEO,ELEM,qsi,w,CDC,k],3)
# H,G=cal_Aeb(Tree[14],Tree[8],[NOS,NOS_GEO,ELEM,qsi,w,CDC,k])

println("A matriz foi dividida em $(length(block[:,3])) blocos")
println("Dentre estes blocos, $(sum(block[:,3])) são aproximados por matrizes de baixo rank")
println("Rank máximo: $(ninterp*ninterp)")
println("norma = $(norm(T-NOS[:,4]))")
