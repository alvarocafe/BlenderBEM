function calc_solfund(x,y,xd,yd,nx,ny,k)
#Calcula as soluções fundamentais
r=sqrt((x-xd)^2+(y-yd)^2); # Raio (distância entre ponto fonte e
                           # ponto campo)
rx=(x-xd); # Componente x do raio
ry=(y-yd); # Componente y do raio
Tast=-1/(2*pi*k)*log(r); # Solução fundamental da temperatura
qast=1/(2*pi)*(rx*nx+ry*ny)/r^2; # Solução fundamental do fluxo
return Tast, qast
end

function calcula_HeGns(x1,y1,x2,y2,xd,yd,qsi,w,k)
#integração não singular
n_pint=length(qsi); # Número de pontos de integração (comprimento do
                  #    vetor qsi)
G=0; # Inicializa o somatorio de g
H=0; # Inicializa o somatorio de h
for kk=1:n_pint # Laço sobre os pontos de integração
    N1,N2=calc_fforma(qsi[kk]); # Calcula as funções de forma
    dN1dqsi= -0.5
    dN2dqsi=0.5 # Calcula as derivadas das
                                                     #    funções de forma
    x=N1*x1+N2*x2; # Calcula a coordenada x do ponto de integração
    y=N1*y1+N2*y2; # Calcula a coordenada y do ponto de integração

    dxdqsi=dN1dqsi*x1+dN2dqsi*x2;
    dydqsi=dN1dqsi*y1+dN2dqsi*y2;
    dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);

    sx=dxdqsi/dgamadqsi; # Componente x do vetor tangente
    sy=dydqsi/dgamadqsi; # Componente y do vetor tangente
    nx=sy; # Componente x do vetor normal
    ny=-sx; # Componente y do vetor normal

    Tast,qast=calc_solfund(x,y,xd,yd,nx,ny,k); # Calcula as soluções fundamentais
    H=H+qast*dgamadqsi*w[kk]; # Integral da matriz H
    G=G+Tast*dgamadqsi*w[kk]; # Integral da matriz G
end
return G, H
end

function calcula_HeGs(x1,y1,x2,y2,k)
# integração singular
  H=-1/2; # Matriz elementar h
  L=sqrt((x2-x1)^2+(y2-y1)^2); # Comprimento do elemento
  G=(L/(2*pi*k))*(1-log(L/2)); # Matriz elementar g
  return G, H
end

function Gauss_Legendre(x1,x2,n)
  eps=3e-14;
  m::Int64 = round((n+1)/2);
  x = zeros(1,n)
  w = zeros(1,n)
  pp = 1
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for i=1:m
    z=cos(pi*(i-0.25)/(n+0.5));
    while 1==1
      p1=1.0;
      p2=0.0;
      for j=1:n
        p3=p2;
        p2=p1;
        p1=((2*j-1)*z*p2-(j-1)*p3)/j;
      end
      pp=n*(z*p1-p2)/(z*z-1);
      z1=z;
      z=z1-p1/pp;
      if(abs(z-z1)<eps)
        break
      end
    end
    x[i]=xm-xl*z;
    x[n+1-i]=xm+xl*z;
    w[i]=2*xl/((1-z*z)*pp*pp);
    w[n+1-i]=w[i];
  end
  return x,w
end
function calc_fforma(qsi)
  # Calcula as funções de forma lineares contínuas N1 e N2
  N1=1. /2. .*(1. .-qsi); # Função de forma N1 => linear contínua
  N2=1. /2. .*(1. .+qsi); # Função de forma N2 => linear contínua
  return N1,N2
end

function monta_Teq(CDC,x)
  # Separa fluxo e temperatura
  # ncdc = n�mero de linhas da matriz CDC
  # T = vetor que cont�m as temperaturas nos n�s
  # q = vetor que cont�m o fluxo nos n�s

  ncdc = length(CDC[:,1]);
  nnos = length(x)
  T = zeros(nnos,1)
  q = zeros(nnos,1)
  for i=1:ncdc # Laco sobre as condicoes de contorno
    tipoCDC=CDC[i,2]; # Tipo da condi��o de contorno
    valorCDC=CDC[i,3]; # Valor da condi��o de contorno
    valorcalculado=x[i]; # Valor que antes era desconhecido
    if tipoCDC == 1 # Fluxo � conhecido
      T[i] = valorcalculado; # A temperatura � o valor calculado
      q[i] = valorCDC; # O fluxo � a condi�ao de contorno
    else # A temperatura � conhecida
      T[i] = valorCDC; # A temperatura � a condi�ao de contorno
      q[i] = valorcalculado; # O fluxo � o valor calculado
    end
  end

  return T,q
end




