function montamatrizvetor(NOS, NOS_GEO, ELEM, k, CDC)
nelem=size(ELEM,1); # N�mero de elementos (n�mero de linhas da
#  matriz ELEM)

nnos=nelem; # N�mero de n�s
G=zeros(nnos,nnos); # Inicializa��o da matriz G
H=zeros(nnos,nnos); # Inicializa��o da matriz H

npg=12; # N�mero de pontos de Gauss
qsi_tri,w_tri=Gauss_Legendre(0,1,npg); # Pontos e pesos de Gauss do
                                         # elemento triangular
qsi,w=Gauss_Legendre(-1,1,npg); # Pontos e pesos de Gauss do elemento
                                  # quadrilateral
A = complex(zeros(nnos,nnos));
b = complex(zeros(nnos, 1));

for j=1:nelem # La�o sobre os elementos
	tipoCDC = CDC[j,2]; #Tipo da condicao de contorno CDC[pc,2] = 0 condicao de pressao, =1 condicao de fluxo
	#valorCDC = complex(CDC[pc,3],CDC[pc,4])	#Valor da condicao no elementos
	valorCDC = CDC[j,3];
	nos = ELEM[j,2:4];		#Nos contidos no elementos

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

        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
                     # normal ao elemento
        for i=1:nnos # La�o sobre os pontos fontes
            xd=NOS[i,2]; # Coordenada x do ponto fonte
            yd=NOS[i,3]; # Coordenada y do ponto fonte
            zd=NOS[i,4]; # Coordenada y do ponto fonte

            if i==j # O ponto fonte pertence ao elemento
                G[i,j],H[i,j]=calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,k); # Integra��o singular
            else # O ponto fonte n�o pertence ao elemento
                G[i,j],H[i,j]=calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi_tri,w_tri,k); # Integra��o regular
            end
        #Aplica as condicoes de contorno
        	if tipoCDC == 0 # A pressao eh conhecida
        		A[i,j] = -G[i,j]; 	#Os valores de G vao para a matriz A
                b[i,1] = b[i,1] - H[i,j]*valorCDC; #Os valores de H vao para o vetor b
            else
                A[i,j] = +H[i,j]; 	#Os valores de H vao para a matriz A
                b[i,1] = b[i,1] + G[i,j]*valorCDC; #Os valores de G vao para o vetor b
            end
        end
end
return A,b
end

