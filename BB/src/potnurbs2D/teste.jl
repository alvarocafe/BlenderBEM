include("decomp.jl")
using Plots,FastGaussQuadrature,LinearAlgebra
# plotly()
npg=10; # Número de pontos que são usados para calcular a curva
#xi=linspace(-1,1,npg) 
xi,wg=gausslegendre(npg)
#xi=[-0.861136 -0.339981 0.339981 0.861136];
#wg=[0.347855 0.652145  0.652145  0.347855];
#npg=length(xi)
knot=[0 0 0 0 1 2 3 4 4 4 4]
p=3 # Grau da B-Spline
P=[0.0  0.0  3.0  7.0  9.0  13.0   16.0
   0.0  1.0  1.0  0.0  2.0   2.25   1.0]
n=size(P,2)-1 # número de pontos de controle menos 1
w=ones(n+1,1) # pesos
# C = operador de extração 
# nb = número de curvas de Bézier
# conn = conectividade da spline (diz quais são os pontos de controle de cada trecho da curva)
C,nb,conn = bezierExtraction(knot,p);
shapes=Vector(npg)
derivs=Vector(npg)
B1=zeros(npg)
B2=zeros(npg)
B3=zeros(npg)
B4=zeros(npg)
kmat=1
# Loop para gerar as funções de forma (polinômios de bernstein) e suas derivadas
for gp=1:npg
    shapes[gp],derivs[gp] = bernsteinbasis(p,0,xi[gp],0)
    B1[gp]=shapes[gp][1]
    B2[gp]=shapes[gp][2]
    B3[gp]=shapes[gp][3]
    B4[gp]=shapes[gp][4]
end
XX=zeros(nb*npg) # Vetor que armazena as coordenadas x da curva
YY=zeros(nb*npg) # Vetor que armazena as coordenadas y da curva
R1=zeros(nb*npg) # Vetor que armazena as coordenadas x da curva
R2=zeros(nb*npg) # Vetor que armazena as coordenadas y da curva
R3=zeros(nb*npg) # Vetor que armazena as coordenadas x da curva
R4=zeros(nb*npg) # Vetor que armazena as coordenadas y da curva
for j=1:nb # Loop sobre as curvas de Bézier (que corresponde aos elementos de contorno)
    for i = 1 : npg # Percorre os pontos de integração
        R,dRdxi=basisfundecomp(shapes[i],derivs[i],C[:,:,j],w[conn[j]])
        xy=P[1:2,conn[j]]*R
        dxydxi=P[1:2,conn[j]]*dRdxi;
        # derivadas das fun��es de forma
        x=xy[1]; # Calcula a coordenada x do ponto de integra��o
        y=xy[2] # Calcula a coordenada y do ponto de integra��o
        XX[npg*(j-1)+i]=x        
        YY[npg*(j-1)+i]=y
        R1[npg*(j-1)+i]=R[1]        
        R2[npg*(j-1)+i]=R[2]
        R3[npg*(j-1)+i]=R[3]        
        R4[npg*(j-1)+i]=R[4]
    end
end
p1=plot(P[1,:],P[2,:],color="blue",label="Control points",marker=(1,0.3,:o))
p1=plot!(XX,YY,label="B-Spline curve",marker=(1,0.3,:x),color="red",legend=:bottomright)


p2=plot(xi,B1,label="B1",marker=(1,0.3,:d),color="red",legend=:topright)
p2=plot!(xi,B2,label="B2",marker=(1,0.3,:o),color="green")
p2=plot!(xi,B3,label="B3",marker=(1,0.3,:s),color="blue")
p2=plot!(xi,B4,label="B4",marker=(1,0.3,:s),color="black")
p2=title!("Polinômios de Berstein de ordem 3")

p3=plot(XX,R1,label="R1",legend=true,marker=(1,0.3,:x),color="red")
p3=plot!(XX,R2,label="R2",legend=true,marker=(1,0.3,:o),color="green")
p3=plot!(XX,R3,label="R3",legend=true,marker=(1,0.3,:s),color="blue")
p3=plot!(XX,R4,label="R4",legend=true,marker=(1,0.3,:s),color="black")
p3=title!("Funções de forma B-Splines racionais")

