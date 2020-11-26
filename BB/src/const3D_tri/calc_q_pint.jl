function calc_q_pint(PONTOS_int,NOS_GEO,ELEM,T,q,k,qsi,w,inc)
# Evaluates the flux of the velocity potential at domain points
n_pint=length(PONTOS_int[:,1]); # Number of domain points
n_elem=length(T); # Number of elements
Gx_int = complex(zeros(n_pint,n_elem))
Gy_int = complex(zeros(n_pint,n_elem))
Gz_int = complex(zeros(n_pint,n_elem))
Hx_int = complex(zeros(n_pint,n_elem))
Hy_int = complex(zeros(n_pint,n_elem))
Hz_int = complex(zeros(n_pint,n_elem))
phi_inc = complex(zeros(n_pint,1))
g = complex(zeros(n_pint,1))
for i=1:n_pint # Loop over the integration points
    x_fonte=PONTOS_int[i,2]; # x coordinate of the source point
    y_fonte=PONTOS_int[i,3]; # y coordinate of the source point
    z_fonte=PONTOS_int[i,4]; # z coordinate of the source point
    for j=1:n_elem  # Loop over the elements
        no1=ELEM[j,2]; # First node of the element
        no2=ELEM[j,3]; # Second node of the element
        no3=ELEM[j,4]; # Third node of the element

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
        Gx_int[i,j],Gy_int[i,j],Gz_int[i,j],Hx_int[i,j],Hy_int[i,j],Hz_int[i,j]=calcula_SGeSHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x_fonte,y_fonte,z_fonte,n,qsi,w,k); # Evaluates the kernel for the hypersingular boundary equation
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
# Velocity potential flux at domain points
dphidx=Hx_int*T-Gx_int*q;
dphidy=Hy_int*T-Gy_int*q; 
dphidz=Hz_int*T-Gz_int*q; 
#T_pint = - (H_int*T - G_int*q - phi_inc)
#T_pint=-(H_int*T'-G_int*q'-g'); 
return dphidx,dphidy,dphidz
end

