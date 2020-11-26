function calc_T_pint(PONTOS_int,NOS_GEO,ELEM,T,q,k,qsi,w,inc)
# Calcula a temperatura nos pontos internos
n_pint=length(PONTOS_int[:,1]); # Numero de pontos internos
n_elem=length(T); # Numero de elementos
G_int = complex(zeros(n_pint,n_elem))
H_int = complex(zeros(n_pint,n_elem))
phi_inc = complex(zeros(n_pint,1))
g = complex(zeros(n_pint,1))
for i=1:n_pint # La�o sobre os pontos internos
    x_fonte=PONTOS_int[i,2]; # Coordenada x do ponto fonte
    y_fonte=PONTOS_int[i,3]; # Coordenada y do ponto fonte
    z_fonte=PONTOS_int[i,4]; # Coordenada z do ponto fonte
    for j=1:n_elem  #La�o sobre os elementos
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
        G_int[i,j],H_int[i,j]=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x_fonte,y_fonte,z_fonte,n,qsi,w,k); # Chama a functio para calculo de H e G
        # quando o ponto fonte nao pertence ao elemento
    end
#    if(inc[1,1]!=0)
#    	#g[i,1]=calc_q(x_fonte,y_fonte,fc,FR,CW,GE);
#			g[i,1]=0;
#		else
#    	g[i,1]=0;
#    end
#		if inc[1,1] != 0
	#Vamos incluir um termo de onda incidente
#				phi_inc[i,1] = calc_inc(x_fonte,y_fonte,z_fonte,FR,CW,inc[1,:]);
#		end
end
T_pint = - (H_int*T - G_int*q - phi_inc)
#T_pint=-(H_int*T'-G_int*q'-g'); # Vetor que contem a temperatura nos
#      pontos internos
return T_pint
end

