#Input dos dados
p1 = [0 0 0]
p2 = [1 1 0]
p3 = [2 0 0]
NOS_GEO,NOS,ELEM,CDC = dad_triconst(p1,p2,p3)


#Geração dos pontos e pesos de Gauss
npg=6;
qsi,w = Gauss_Legendre(0,1,npg)

#Integração dos elementos triangulares constantes
I1 = complex(0,0)
I2 = complex(0,0)
x = zeros(size(qsi,2),size(qsi,2))
y = zeros(size(qsi,2),size(qsi,2))
z = zeros(size(qsi,2),size(qsi,2))
N = zeros(size(qsi,2),3)
dN = zeros(size(qsi,2),6)
dgamadqsi = zeros(size(qsi,2),size(qsi,2))
phiast =complex(zeros(size(qsi,2),size(qsi,2)))
qast = complex(zeros(size(qsi,2),size(qsi,2)))
for i = 1:size(qsi,2)
  for j = 1:size(qsi,2)
    xi = (1 - qsi[j])*qsi[i]
    N[j,:] =calc_fforma_triconst(xi,qsi[j]); # Calcula as funções de forma
    dN[j,:] = calc_dfforma_triconst(xi,qsi[j])
    x[i,j]=N[j,1]*NOS_GEO[1,2]+N[j,2]*NOS_GEO[2,2]+N[j,3]*NOS_GEO[3,2]; # Calcula a coordenada x do ponto de integração
    y[i,j]=N[j,1]*NOS_GEO[1,3]+N[j,2]*NOS_GEO[2,3]+N[j,3]*NOS_GEO[3,3]; # Calcula a coordenada y do ponto de integração
    z[i,j]=N[j,1]*NOS_GEO[1,4]+N[j,2]*NOS_GEO[2,4]+N[j,3]*NOS_GEO[3,4]; # Calcula a coordenada y do ponto de integração
    dgamadqsi[i,j]=cal_Jacobiano3D(NOS_GEO[:,2], NOS_GEO[:,3], [0 0 0], dN[j,:]); # Jacobiano
    #dx=dN[kk,1]*x[1]+dN[kk,2]*x[2]+dN[kk,3]*x[3]; # Calcula a coordenada x do ponto de integração
    #dy=dN[kk,1]*y[1]+dN[kk,2]*y[2]+dN[kk,3]*y[3]; # Calcula a coordenada y do ponto de integração
    #sx=dx/dgamadqsi[kk]; # Componente x do vetor tangente
    #sy=dy/dgamadqsi[kk]; # Componente y do vetor tangente
    #nx=sy; # Componente x do vetor normal
    #ny=-sx; # Componente y do vetor normal
    #phiast[kk],qast[kk] =calc_solfund(xx[kk],yy[kk],x[1],y[1],nx,ny,340,340); # Obtemos as soluções fundamentais
    I1=I1+1*w[i]*w[j]*dgamadqsi[i,j]; # Integral da matriz H
    I2=I2+1*w[i]*w[j]*dgamadqsi[i,j]; # Integral da matriz G
  end
end
println("O valor das integrais I1 e I2 são: I1= ",I1," e I2= ",I2)

#Iniciamos a parte gráfica
close("all") #Vamos fechar as janelas que possam estar abertas

#Solução Fundamental
#figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
#plot(qsi',abs(phiast),label = "phiast",linestyle="-",linewidth=3,color ="green")
#plot(qsi',abs(qast),label = "qast",linestyle="-",linewidth=3,color="blue")
#legend(fontsize="12.0",loc="best")
#xlabel("xi")
#ylabel("Fundamental solutions")
#grid(1)
#title("Fundamental solution for the linear element")
#savefig("docs/discretização/figuras/fundamental_const.pdf",transparent = true)


#Funções de forma
figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
plot(qsi',N[:,1],label = "N1",linestyle="-",linewidth=3,color ="green")
plot(qsi',N[:,2],label = "N2",linestyle="-",linewidth=3,color ="blue")
plot(qsi',N[:,3],label = "N3",linestyle="-",linewidth=3,color ="red")
legend(fontsize="12.0",loc="best")
xlabel("xi")
ylabel("Shape functions")
grid(1)
title("Shape functions for the triangular linear element")
#savefig("docs/discretização/figuras/ff_forma_triconst.pdf",transparent = true)

#Jacobiano
#figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
#plot(qsi',dgamadqsi,label = "Jacobian",linestyle="-",linewidth=3,color ="red")
#surf(qsi',qsi',dgamadqsi,label = "Jacobian")
#legend(fontsize="12.0",loc="best")
#xlabel("xi")
#ylabel("eta")
#zlabel("Jacobian")
#grid(1)
#title("Jacobian for the triangular constant element")
#savefig("docs/discretização/figuras/jacobian_const.pdf",transparent = true)


#Elemento
figure(figsize=[8,4]) #Abrimos uma nova figura com dimensões 8x4 in
plot(NOS_GEO[:,2],NOS_GEO[:,3],linestyle="none",linewidth=3,color="grey",marker="o",markersize=12,label="Geometrical points")
plot(x,y,marker=".",color="black",linestyle="--",markersize=12,label = "Gauss' points")
plot(NOS_GEO[:,2],NOS_GEO[:,3],marker="s",markersize=10,linestyle="none",color="blue",label = "Physical points (Nodes)")
grid(1)
xlabel("x",fontsize="12.0")
ylabel("y",fontsize="12.0")
title("Triangular linear element",fontsize="16.0")
#legend(fontsize="14.0",loc="best")
PyPlot.xlim(-0.2, 2.2)
PyPlot.ylim(-0.2, 1.2)
#savefig("docs/discretização/figuras/triconstante.pdf",transparent = true)
