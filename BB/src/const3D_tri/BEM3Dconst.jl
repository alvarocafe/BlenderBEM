# Boundary element method implementation for the Helmholtz equation using constant bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Modules necessary: SpecialFunctions.jl
using SpecialFunctions
using PyPlot
include("dep.jl") # Includes the dependencies
include("dad4.jl")
FR = 774*100 # Frequency of the problem [Hz]
CW = 343*1000 # Wave propagation speed [mm/s]
k = FR/CW # Wave number
# Gaussian quadrature - generation of points and weights [-1,1]
npg=4; # Number of integration points
qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
println("Importing mesh...")
#@time NOS_GEO,ELEM,elemint,CDC = lermsh("../../dados/vocal_tract.msh",3) #Read the mesh generated using Gmsh
#CCFace = ones(28,3);
#for i = 1:28
#	CCFace[i,:] = [i 1 1]
#end
NOS_GEO, ELEM, CCFace, fc, k = dad4()
#FR = k
#CW = k
#NOS_GEO = [1 0 0 0
#	   2 1 0 0
#	   3 1 1 0
#	   4 0 1 0]
#ELEM = [1 1 2 3 1
#	2 2 3 4 1]
NOS = mostra_geoTRI(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
CDC = gera_CDC(ELEM,CCFace); #Monta a matriz de condicoes de contorno
n_pint = 100
L = 1;
PONTOS_int = zeros(n_pint,4)
for i = 1:n_pint
	PONTOS_int[i,:] = [i 0 0 (L/n_pint)*i]
end
println("Building G and H matrices...")
@time G,H,phi_inc = cal_GeH_POT(NOS,NOS_GEO,ELEM,k,qsi,w,0) #Compute the G and H matrices
println("Applying boundary conditions to build A and b for the linear system...")
@time A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
println("Solving the linear system...")
@time x = A\b # Solves the linear system
println("Separating acoustic pressure from flux...")
@time phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
println("Solving for domain points.")
T_pint=calc_T_pint_POT(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,fc);
q_pint=calc_q_pint_POT(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,fc);

