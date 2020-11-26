function calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k)
#integra��o n�o singular
n_pint=length(qsi); # N�mero de pontos de integra��o.
g=complex(0,0); # Inicializa o somatorio de g
h=complex(0,0); # Inicializa o somatorio de h
N = zeros(3)
for l=1:n_pint # La�o sobre os pontos de integra��o
    for m=1:n_pint # La�o sobre os pontos de integra��o
	xi = (1 - qsi[l])*qsi[m]
        N =calc_fformatri(xi,qsi[l]); #  fun��es de forma
	#dN = calc_dfforma(xi,qsi[m])
        x=N[1]*x1+N[2]*x2+N[3]*x3; # coordenada x do ponto de integra��o
        y=N[1]*y1+N[2]*y2+N[3]*y3 # coordenada y do ponto de integra��o
        z=N[1]*z1+N[2]*z2+N[3]*z3; # coordenada z do ponto de integra��o
	J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,xi,qsi[m]);# jacobiano
        #Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,FR,CW); # Solu��es
        #  fundamentais
	Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,k); # Solu��es
        #  fundamentais
	#Tast = 1 #To test the numerical integration, the value of the fundamental solutions is set to 1. In this special case, the value of the sum of the rows of matrices H or G should equal the surface area of the model.
	#qast = 1        
	g=g+Tast*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz G
        h=h+qast*complex((1-qsi[l])*w[l]*w[m]*J,0); # Integral da matriz H
    end
end
return g,h
end

