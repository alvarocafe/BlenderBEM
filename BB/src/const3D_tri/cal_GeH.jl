function cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,inc)
# Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.

nelem=size(ELEM,1); # N�mero de elementos (n�mero de linhas da
#  matriz ELEM)
qsitelles,Jtelles = telles(qsi,0); # Evaluates the Telles' points and Jacobian
npg=12; # Number of integration points
qsiquad,wquad = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
nnos=nelem; # N�mero de n�s
G=complex(zeros(nnos,nnos)); # Inicializa��o da matriz G
H=complex(zeros(nnos,nnos)); # Inicializa��o da matriz H
phi_inc = complex(zeros(nelem,1));
for i=1:nnos # La�o sobre os pontos fontes
		xd=NOS[i,2]; # Coordenada x do ponto fonte
		yd=NOS[i,3]; # Coordenada y do ponto fonte
		zd=NOS[i,4]; # Coordenada y do ponto fonte

		for j=1:nelem # La�o sobre os elementos
		    no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
		    no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
		    no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

		    x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
		    y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
		    z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

		    x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
		    y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
		    z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

		    x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
		    y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
		    z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

		    n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio normal ao elemento
		        if i==j # O ponto fonte pertence ao elemento
		           #G[i,j],H[i,j]=calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,k); # Integra��o singular
			   G[i,j],H[i,j]= calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,k)# Integra��o singular
			   #G[i,j],H[i,j]=calcula_GeHns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o singular
			    #Gtelles,Htelles=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsitelles,w.*Jtelles,k); # Integra��o singular
#G[i,j]=1
#H[i,j]= -0.5
#			erro = (G[i,j] - Gtelles)
#			println("diferença entre g e gtelles= $erro")
		        else # O ponto fonte n�o pertence ao elemento
		           # G[i,j],H[i,j]=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o
		            #  regular
			   G[i,j],H[i,j]=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o singular
		        end
		    end
				if inc[1,1] != 0
			#Vamos incluir um termo de onda incidente
					phi_inc[i,1] = calc_inc(xd,yd,zd,k,k,inc[1,:]);
				end
end

return G,H,phi_inc
end

