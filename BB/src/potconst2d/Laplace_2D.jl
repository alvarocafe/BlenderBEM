include("dad.jl")
include("hmat.jl")
include("MEC.jl")
include("pre_proc.jl")
include("tree.jl")
include("inpoly.jl")
include("mostra_problema.jl")
using Statistics, LinearAlgebra, KrylovMethods

PyPlot.close("all")   # close all plot windows

PONTOS, SEGMENTOS, MALHA, CCSeg, k = dad_0()
NOS_GEO,NOS,ELEM,CDC = format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg)

nnos=size(NOS,1);
println("Número de nós = $nnos")
npg=8;
qsi,w = Gauss_Legendre(-1,1,npg) # Gera os pontos e pesos de Gauss

b1 = 1:1:nnos
Ac,Bc = cal_Aeb(b1,b1,[NOS,NOS_GEO,ELEM,qsi,w,CDC,k])
bc = Bc*CDC[:,3]
println(sum(sum(bc)))
xc = bc\Ac
Tc,qc = monta_Teq(CDC,xc)

Tree,child,center_row,diam,inode,ileaf = cluster(NOS[:,2:3],floor(sqrt(2*length(NOS[:,1]))))

ninterp=4 # Número de pontos de interpolação
η =.4 # Coeficiente relacionado a admissibilidade dos blocos

allow=checa_admiss(η,center_row,diam,inode,ileaf) # mostra quais blocos são admissíveis

block = blocks(Tree,child,allow) # Função que retorna os blocos admissiveis


hmati,bi =Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,qsi,w,CDC,k],ninterp)

xi = gmres(vet->matvec(hmati,vet,block,Tree),bi,5,tol=1e-8,maxIter=1000,out=0) 

T,q=monta_Teq(CDC,xi[1])

# A=montacheia(hmati,block,Tree,nnos)
# L,HH,G2=cal_Aeb_interp(Tree[14],Tree[8],[NOS,NOS_GEO,ELEM,qsi,w,CDC,k],3)
# H,G=cal_Aeb(Tree[14],Tree[8],[NOS,NOS_GEO,ELEM,qsi,w,CDC,k])

println("A matriz foi dividida em $(length(block[:,3])) blocos")
println("Dentre estes blocos, $(sum(block[:,3])) são aproximados por matrizes de baixo rank")
println("Rank máximo: $(ninterp*ninterp)")


mostra_problema(ELEM,NOS_GEO,NOS)


[NOS[:,2] T Tc]
