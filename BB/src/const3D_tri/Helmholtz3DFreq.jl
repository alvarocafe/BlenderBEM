# Programa de elementos de contorno aplicado acustica tri-dimensioanal no
# dominio da frequencia
# Tipo de elementos: constantes
# Por: Eder Lima de Albuquerque e Alvaro Campos Ferreira
# Brasilia, dezembro de 2016
workspace()
using PyPlot
using PyCall
@pyimport matplotlib.colors as col
@pyimport matplotlib.cm as cm
@pyimport mpl_toolkits.mplot3d as mp
@pyimport mpl_toolkits.mplot3d.art3d as ar
plt=PyPlot

#include("dad3_pot.jl") # Arquivo de entrada de dados
include("BEM3D_fix.jl")
include("lermsh.jl")	# Leitura de arquivos de malha criados pelo GMSH

NOS_GEO,ELEM,elemint,CDC = lermsh("../../dados/cilindro_10x100mm.msh",3)

# Mostra a geometria do problema
#figure(1)
CCFace = [1 1 10
	  2 0 10
	  3 0 10]
CDC = gera_CDC(ELEM,CCFace); #Monta a matriz de condicoes de contorno
NOS = mostra_geo(NOS_GEO,ELEM);  #Monta a matriz que contem os nos fisicos
#NFR = length(AFR); # Numero de frequencias que serao analisadas
NFR = 1
AFR = 10000
CW = 343*1000
k=AFR/CW
npint = 100
PONTOS_int = zeros(npint,4)
for i = 1:npint
	PONTOS_int[i,:] = [i 0.0 0.0 (100/npint)*i]
end
println("2. Resolvendo o problema para ", (NFR)," frequencias");
nnos=size(NOS,1);
pfreq=complex(zeros(nnos,NFR));
qfreq=complex(zeros(nnos,NFR));
pint_freq=complex(zeros(npint,NFR));

npg=12; # N�mero de pontos de Gauss
qsi,w = Gauss_Legendre(0,1,npg); # Pontos e pesos de Gauss
npg=12; # N�mero de pontos de Gauss
qsi_quad,w_quad = Gauss_Legendre(-1,1,npg); # Pontos e pesos de Gauss


#Onda incidente
inc = [1 1 0 0 0.1];
A = complex.(zeros(nnos,nnos))
for  i=1:NFR
    println("Frequencia ", i ," = ", AFR[i] ," rad/s");
    FR=AFR[i];
#    G, H,phi_inc = cal_HeG(NOS, NOS_GEO, ELEM,FR,CW,qsi,w,inc); # Monta as matrizes H e G
#    A,b= aplica_cdc(G,H,CDC); # Aplica as condicoes de contorno
    A,b =  monta_matrizvetor(NOS,NOS_GEO, ELEM, FR, CW, qsi, w, qsi_quad,w_quad, inc);
    A1,b1 = montamatrizvetor(NOS,NOS_GEO,ELEM,k,CDC);
#    x = A\(b+phi_inc); # Calcula as variaveis desconhecidas
    error = norm(A-A1);
    println(error)
#    println(A1)
    x = A\(b); # Calcula as variaveis desconhecidas
    p,q= monta_Teq(CDC,x); # Monta T e q
    p_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,p,q,CW,FR,qsi,w,inc);
    pfreq[:,i]=(p);
    qfreq[:,i]=(q);
    #pint_freq[:,i]=(p_pint);
end
mostra_resultados2(NOS_GEO,ELEM,real(pfreq[:,1]))
#println("Pressão acústica= ",pfreq[:,1])
