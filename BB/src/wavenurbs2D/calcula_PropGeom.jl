function calcula_PropGeom(collocCoord,nnos,crv,dcrv,kmat)
dcrvs=map(x->nrbderiv(x),crv)
n = length(crv);	# Number of curves
ncollocpoints=size(collocCoord,1)

H=zeros(ncollocpoints,ncollocpoints);
G=zeros(ncollocpoints,ncollocpoints);
npgauss=12;
qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
for i = 1 : n

uk=unique(crv[i].knots)

ne=length(uk)-1;
for j=1:ne
  range=uk[j+1]-uk[j]
  mid=(uk[j+1]+uk[j])/2
  for k=1:ncollocpoints
          xfonte=collocCoord[k,1]
          yfonte=collocCoord[k,2]


    # Integrais de dom�nio
    g,h,id=integra_elem(xfonte,yfonte,crv[i],dcrv[i],range,mid,qsi,w,kmat); # Integra��o sobre o
                               #  element (I = int F n.r/r dGama)
    H[k,id+nnos[i]]=H[k,id+nnos[i]]+h;
    G[k,id+nnos[i]]=G[k,id+nnos[i]]+g;
end
end
end
return H,G
end
function  integra_elem(xfonte,yfonte,crv,dcrv,range,mid,qsi,w,k)
C,nb,conn=bezierExtraction(crv.knots,crv.order-1)

# Integra��o sobre os elementos (integral I)

# N�mero de pontos de Gauss usados na integra��o de I
npg = length(qsi);
shapes=Vector(npg)
derivs=Vector(npg)
for gp=1:npg
        shapes[gp],derivs[gp] = bernsteinbasis(crv.order-1,0,qsi[gp],0)
end


dudqsi=range/2
g = zeros(crv.order)
h = zeros(crv.order)
id=zeros(Int,crv.order,1)
for i = 1 : npg # Percorre os pontos de integra��o
    p,dp=basisfundecomp(shapes[i],derivs[i],C[j],crv.coefs[3,conn[j]],crv.coefs[1:2,conn[j]])
       # derivadas das fun��es de forma
    dgamadu=norm(dp)
    x=p[1]-xfonte # Calcula a coordenada x do ponto de integra��o
    y=p[2]-yfonte # Calcula a coordenada y do ponto de integra��o
    dxdqsi=dp[1];
    dydqsi=dp[2];
    dgamadqsi=√(dxdqsi^2+dydqsi^2);

    sx=dxdqsi/dgamadqsi; # Component x do vetor tangente
    sy=dydqsi/dgamadqsi; # Componente y do vetor tangente

    nx=sy; # Componente x do vetor normal unit�rio
    ny=-sx; # Componente y do vetor normal unit�rio



    r=√(x^2+y^2); # raio
    rx=x/r; # Componente x do vetor unit�rio r
    ry=y/r; # Componente y do vetor unit�rio r
    nr=nx*rx+ny*ry; # Produto escalar dos vetores unit�rios r e n
    Tast=-1/(2*pi*k)*log(r)
    qast=1/(2*pi)*(rx*nx+ry*ny)/r
    h = h + shapes[i]*qast*dgamadu *dudqsi* w[i]
    g = g + shapes[i]*Tast*dgamadu *dudqsi* w[i]
end
return g,h,id
end

function aplica_CDC(G,H,CDC,E)
# Aplica as condições de contorno trocando as colunas das matrizes H e G

ncdc = size(CDC,1); # número de linhas da matriz CDC
A=H;
B=G;
for i=1:ncdc # Laço sobre as condições de contorno
    tipoCDC = CDC[i,2]; # Tipo da condição de contorno
    if tipoCDC == 0 # A temperatura é conhecida
        colunaA=-A[:,i]; # Coluna da matriz H que será trocada
        A[:,i]=-B[:,i]; # A matriz H recebe a coluna da matriz G
        B[:,i]=colunaA; # A mstriz G recebe a coluna da matriz H
    end
end
valoresconhecidos=E\CDC[:,3] # Valores das condições de contorno
b=B*valoresconhecidos; # vetor b
return A,b
end

function monta_Teq(CDC,x)
# Separa fluxo e temperatura

# ncdc = número de linhas da matriz CDC
# T = vetor que contém as temperaturas nos nós
# q = vetor que contém o fluxo nos nós

ncdc = size(CDC,1);
T=zeros(ncdc)
q=zeros(ncdc)
for i=1:ncdc # Laço sobre as condições de contorno
    tipoCDC=CDC[i,2] # Tipo da condição de contorno
    valorCDC=CDC[i,3] # Valor da condição de contorno
    valorcalculado=x[i] # Valor que antes era desconhecido
    if tipoCDC == 1 # Fluxo é conhecido
        T[i] = valorcalculado; # A temperatura é o valor calculado
        q[i] = valorCDC; # O fluxo é a condiçao de contorno
    else # A temperatura é conhecida
        T[i] = valorCDC; # A temperatura é a condiçao de contorno
        q[i] = valorcalculado; # O fluxo é o valor calculado
    end
end
return T,q
end
