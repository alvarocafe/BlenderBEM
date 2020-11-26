function calcula_Gs(X1,X2,X3,Xd,W,FR,CW,N1,N2,N3,N4,dN1dqsi,dN2dqsi,dN3dqsi,dN4dqsi,dN1deta,dN2deta,dN3deta,dN4deta)
	#Essa eh a funcao que calcula a integral singular da matriz G
	g = 0.0; # inicialização da matriz G

	for kk=1:3
	    x1t=Xd[1]; # coordenada x do primeiro nó do quadrilatero desgenerado
	    y1t=Xd[2]; # coordenada y do primeiro nó do quadrilatero desgenerado
	    z1t=Xd[3]; # coordenada z do primeiro nó do quadrilatero desgenerado
	    x2t=Xd[1]; # coordenada x do segundo nó do quadrilatero desgenerado
	    y2t=Xd[2]; # coordenada y do segundo nó do quadrilatero desgenerado
	    z2t=Xd[3]; # coordenada z do segundo nó do quadrilatero desgenerado

	    if(kk==1) # Terceiro e quarto nós do primeiro quadrilátero desgenerado
	        x3t=X1[1]; # coordenada x do terceiro nó do quadrilatero desgenerado
	        y3t=X1[2]; # coordenada y do terceiro nó do quadrilatero desgenerado
	        z3t=X1[3]; # coordenada z do terceiro nó do quadrilatero desgenerado
	        x4t=X2[1]; # coordenada x do quarto nó do quadrilatero desgenerado
	        y4t=X2[2]; # coordenada y do quarto nó do quadrilatero desgenerado
	        z4t=X2[3]; # coordenada z do quarto nó do quadrilatero desgenerado

	    elseif(kk==2) # Terceiro e quarto nós do segundo quadrilátero
	                                  # desgenerado
	        x3t=X2[1]; # coordenada x do terceiro nó do quadrilatero desgenerado
	        y3t=X2[2]; # coordenada y do terceiro nó do quadrilatero desgenerado
	        z3t=X2[3]; # coordenada z do terceiro nó do quadrilatero desgenerado
	        x4t=X3[1]; # coordenada x do quarto nó do quadrilatero desgenerado
	        y4t=X3[2]; # coordenada y do quarto nó do quadrilatero desgenerado
	        z4t=X3[3]; # coordenada z do quarto nó do quadrilatero desgenerado

	    elseif(kk==3) # Terceiro e quarto nós do terceiro quadrilátero
	                                  # desgenerado
	        x3t=X3[1]; # coordenada x do terceiro nó do quadrilatero desgenerado
	        y3t=X3[2]; # coordenada y do terceiro nó do quadrilatero desgenerado
	        z3t=X3[3]; # coordenada z do terceiro nó do quadrilatero desgenerado
	        x4t=X1[1]; # coordenada x do quarto nó do quadrilatero desgenerado
	        y4t=X1[2]; # coordenada y do quarto nó do quadrilatero desgenerado
	        z4t=X1[3]; # coordenada z do quarto nó do quadrilatero desgenerado

	    end

	    X1t = [x1t, y1t, z1t];
	    X2t = [x2t, y2t, z2t];
	    X3t = [x3t, y3t, z3t];
	    X4t = [x4t, y4t, z4t];

	    x = N1*X1t[1] + N2*X2t[1] + N3*X3t[1] + N4*X4t[1];
	    y = N1*X1t[2] + N2*X2t[2] + N3*X3t[2] + N4*X4t[2];
	    z = N1*X1t[3] + N2*X2t[3] + N3*X3t[3] + N4*X4t[3];


	    # Jacobiano
			dxdqsi = X1t[1]*dN1dqsi + X2t[1]*dN2dqsi + X3t[1]*dN3dqsi + X4t[1]*dN4dqsi;
	    dydqsi = X1t[2]*dN1dqsi + X2t[2]*dN2dqsi + X3t[2]*dN3dqsi + X4t[2]*dN4dqsi;
	    dzdqsi = X1t[3]*dN1dqsi + X2t[3]*dN2dqsi + X3t[3]*dN3dqsi + X4t[3]*dN4dqsi;

	    dxdeta = X1t[1]*dN1deta + X2t[1]*dN2deta + X3t[1]*dN3deta + X4t[1]*dN4deta;
	    dydeta = X1t[2]*dN1deta + X2t[2]*dN2deta + X3t[2]*dN3deta + X4t[2]*dN4deta;
	    dzdeta = X1t[3]*dN1deta + X2t[3]*dN2deta + X3t[3]*dN3deta + X4t[3]*dN4deta;

	    g1 = dydqsi.*dzdeta - dzdqsi.*dydeta;
	    g2 = dzdqsi.*dxdeta - dxdqsi.*dzdeta;
	    g3 = dxdqsi.*dydeta - dydqsi.*dxdeta;
	    J = real(sqrt.(g1.^2.0 + g2.^2.0 + g3.^2.0));

	    # Solução fundamental

	    r1 = x - Xd[1];
	    r2 = y - Xd[2];
	    r3 = z - Xd[3];
	    R = real(sqrt.(r1.^2 + r2.^2 + r3.^2));
#			ZW=complex(0.,-FR*R/CW);
#			past=1; #exp(ZW)/complex(4*pi*R,0.); #Solucao fundamental para a pressao acustica
			Tast = 1.0./(4.0*FR*pi*R); #Solucao fundamental para a temperatura
			g = g + sum(sum(Tast.*W.*J));

	end
return g
end

