function calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsil,w,k)
#integra��o n�o singular
n_pint=length(qsil); # N�mero de pontos de integra��o.
g=0; # Inicializa o somatorio de g
h=0; # Inicializa o somatorio de h

eta = qsil;
rho = w;
Tast = complex(0,0);
qast = complex(0,0);

for l=1:n_pint # La�o sobre os pontos de integra��o
    for m=1:n_pint # La�o sobre os pontos de integra��o

        qsi=(1-eta[l])*qsil[m];

        N1,N2,N3=calc_fformatri(qsi,eta[l]); #  fun��es de forma
        x=N1*x1+N2*x2+N3*x3; # coordenada x do ponto de integra��o
        y=N1*y1+N2*y2+N3*y3; # coordenada y do ponto de integra��o
        z=N1*z1+N2*z2+N3*z3; # coordenada z do ponto de integra��o

        dNdqsi = [1; 0; -1]; # Derivadas das fun��es de forma
        dNdeta = [0; 1; -1];

        dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3];
        dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3];
        dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3];

        dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3];
        dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3];
        dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3];

        g1 = dydqsi*dzdeta - dzdqsi*dydeta;
        g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
        g3 = dxdqsi*dydeta - dydqsi*dxdeta;
        J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);

        Tast,qast=calc_solfund_POT(x,y,z,xd,yd,zd,n,k); # Solu��es fundamentais

        h=h+qast*(1-eta[l])*rho[l]*w[m]*J; # Integral da matriz H
        g=g+Tast*(1-eta[l])*rho[l]*w[m]*J; # Integral da matriz G

    end
end
return g,h
end

