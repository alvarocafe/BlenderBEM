function monta_matrizvetor(NOS,NOS_GEO, ELEM, FR, CW, qsi, w, qsi_quad,w_quad, inc)
nelem=size(ELEM,1); #Numero de elementos
nnos = nelem; #No caso de elementos constantes, o numero de nos eh o mesmo que o de elementos
A = complex.(zeros(nnos,nnos));
b = complex.(zeros(nnos, 1));
#Inicia os pontos e pesos de Gauss para os elementos triangulares lineares
W = w'*w;
eta = qsi;
ETA = ones(size(qsi))'*qsi;
#Calcula as funcoes de forma triangulares
N1_tri = ((1-eta)'*qsi)';
N2_tri = ones(length(w),1)*eta;
N3_tri = 1 - N1_tri - N2_tri;
#Para a integracao singular, vamos desgenerar em um elemento quadrilateral
W_quad = w_quad'*w_quad;
#Calcula as funcoes de forma quadrilaterais
N1 = (1. /4.)*(1. - qsi_quad')*(1. - qsi_quad);
N2=(1. /4.)*(1. + qsi_quad')*(1. - qsi_quad);
N3=(1. /4.)*(1. + qsi_quad')*(1. + qsi_quad);
N4=(1. /4.)*(1. - qsi_quad')*(1. + qsi_quad);

dN1dqsi = (1. /4.)*ones(npg,1)*(-(1. - qsi_quad));
dN2dqsi = (1. /4.)*ones(npg,1)*(1. - qsi_quad);
dN3dqsi = (1. /4.)*ones(npg,1)*(1. + qsi_quad);
dN4dqsi = (1. /4.)*ones(npg,1)*(-(1. + qsi_quad));

dN1deta = (1. /4.)*ones(npg,1)*(-(1. - qsi_quad));
dN2deta = (1. /4.)*ones(npg,1)*(-(1. + qsi_quad));
dN3deta = (1. /4.)*ones(npg,1)*(1. + qsi_quad);
dN4deta = (1. /4.)*ones(npg,1)*(1. - qsi_quad);

# Inicio da integracao

for pc = 1:nelem #Laco sobre os elementos
#println("pc= ",pc)
	tipoCDC = CDC[pc,2] #Tipo da condicao de contorno CDC[pc,2] = 0 condicao de pressao, =1 condicao de fluxo
	#valorCDC = complex(CDC[pc,3],CDC[pc,4])	#Valor da condicao no elementos
	valorCDC = CDC[pc,3];
	nos = ELEM[pc,2:4]		#Nos contidos no elementos

	X1 = NOS_GEO[nos[1],2:4] 	#Coordenada do primeiro no
	X2 = NOS_GEO[nos[2],2:4] 	#Coordenada do primeiro no
	X3 = NOS_GEO[nos[3],2:4] 	#Coordenada do primeiro no
	v1 = X3 - X2; 	#vetor formado pela aresta 3-2 do elemento
	v2 = X1 - X2; #vetor formado pela aresta 1-2 do elemento
	n = [v1[2]*v2[3] - v2[2]*v1[3]; v2[1]*v1[3]-v1[1]*v2[3]; v1[1]*v2[2] - v2[1]*v1[2]];
	n = n./norm(n);
	J =     J=real(sqrt((-(X1[2]*X2[1])+X1[1]*X2[2]+X1[2]*X3[1]-X2[2]*X3[1] - X1[1]*X3[2]+X2[1]*X3[2])^2+(X1[3]*X2[1]-X1[1]*X2[3]-X1[3]*X3[1]+ X2[3]*X3[1]+X1[1]*X3[3]-X2[1]*X3[3])^2+(-(X1[3]*X2[2])+ X1[2]*X2[3]+X1[3]*X3[2]-X2[3]*X3[2]-X1[2]*X3[3]+X2[2]*X3[3])^2));
	x = N1_tri*X1[1] + N2_tri*X2[1] + N3_tri*X3[1]; #Coordenadas x para o elemento
	y = N1_tri*X1[2] + N2_tri*X2[2] + N3_tri*X3[2]; #Coordenadas x para o elemento
	z = N1_tri*X1[3] + N2_tri*X2[3] + N3_tri*X3[3]; #Coordenadas x para o elemento

	g = Array{Any}
	h = Array{Any}
	for pf = 1: nnos 	#Laco sobre os nos
#println("pf= ",pf)
		Xd = NOS[pf,2:4]; #Coordenadas do ponto fonte
		if pf == pc # Para elementos constantes, se o no corresponder ao elemento, a integracao eh singular
			g = calcula_Gs(X1,X2,X3,Xd,W_quad,FR,CW,N1,N2,N3,N4,dN1dqsi,dN2dqsi,dN3dqsi,dN4dqsi,dN1deta,dN2deta,dN3deta,dN4deta);
			h = -1/2;
		#println("singular")
		else
		#println("não singular")
			r1 = x - Xd[1];
			r2 = y - Xd[2];
			r3 = z - Xd[3];
			R = real(sqrt.(r1.^2 + r2.^2 + r3.^2));
			ZW=complex.(0.,-FR*R/CW);
#			past=1; #exp(ZW)./complex(4*pi*R,0.);
#			drdn=(r1*n[1] + r2*n[2] + r3*n[3])./R;
#			qast=1; #(ZW-complex(1.,0.))*past*complex(drdn./R,0.);	#Solucao fundamental para o fluxo da pressao acustica
			Tast = 1.0./(4.0*FR*pi*R); #Solucao fundamental para a temperatura
			qast = (r1*n[1] + r2*n[2] + r3*n[3])./(4.0*pi*R.^3); #Solucao fundamental para o fluxo da temperatura
			g = sum(sum(Tast.*(1-ETA).*W.*J));
			h = sum(sum(qast.*(1-ETA).*W.*J));
		end
	if tipoCDC == 0 # A pressao eh conhecida
		A[pf,pc] = -g; 	#Os valores de G vao para a matriz A
		b[pf,1] = b[pf,1] - h*valorCDC; #Os valores de H vao para o vetor b
#println("CDC = 0. A[",pf,",",pc,"]= ",A[pf,pc])
#println("CDC = 0. b[",pf,",1]= ",b[pf,1])
	else
			A[pf,pc] = +h; 	#Os valores de H vao para a matriz A
			b[pf,1] = b[pf,1] + g*valorCDC; #Os valores de G vao para o vetor b
	#		println("CDC = 1. A[",pf,",",pc,"]= ",A[pf,pc])
		#	println("CDC = 1. b[",pf,",1]= ",b[pf,1])
	end
end
end
return A,b
end

