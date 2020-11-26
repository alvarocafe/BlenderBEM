function calcula_arco(x1,y1,x2,y2,xc,yc,raio)
# Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
# horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
# and the horizontal direction


dx1 = x1 - xc; dy1 = y1 - yc
dx2 = x2 - xc; dy2 = y2 - yc
tet1=atan(dy1,dx1);
a=[dx1,dy1,0];
b=[dx2,dy2,0];
angle = atan(norm(cross(a,b)),dot(a,b));
if(raio>0)
    tet2=tet1+angle;
else
    tet2=tet1-angle;
end
[tet1,tet2]
end


function calcula_centro(x1,y1,x2,y2,raio)
# Compute the center of an arc given two points and the radius

xm=(x1+x2)/2
ym=(y1+y2)/2
b=√((x2-x1)^2+(y2-y1)^2)
t1=(x2-x1)/b
t2=(y2-y1)/b
n1=t2
n2=-t1
h=√(abs(raio^2-(b/2)^2))
if(raio>0)
   if(n1==0)
      xc=xm
      yc=ym-n2/abs(n2)*h
   else
      xc=-n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
      yc=n2/n1*(xc-xm)+ym
   end
else
   if(n1==0)
      xc=xm
      yc=ym+n2/abs(n2)*h
   else
      xc=n1/abs(n1)*√(h^2*n1^2/(n1^2+n2^2))+xm;
      yc=n2/n1*(xc-xm)+ym
   end
end
[xc,yc]
end

function calc_ncont(SEGMENTOS)
num_seg=size(SEGMENTOS,1);
t=1;
p2=0;
p_ini = SEGMENTOS[1,2];
icont=1; # índice de contornos
contorno=zeros(Int32,1,2)
contorno[1,1]=p_ini;
while(t<num_seg)    # While over all lines
    while(p2!=p_ini)
        p1  = SEGMENTOS[t,2]
        p2  = SEGMENTOS[t,3]
        if p2 == p_ini
            if t < num_seg
                p_ini = SEGMENTOS[t+1,2]
                icont=icont+1;
                contorno=[contorno; zeros(Int32,1,2)]
                contorno[icont,1]=p_ini;
                contorno[icont-1,2]=p_ini-contorno[icont-1,1]
            end;
        end;
        t=t+1
    end                             #end of while p2
end
if(icont>1)
    contorno[icont,2]=num_seg-sum(contorno[1:icont-1,2]);
else
    contorno[1,2]=num_seg;
end
return contorno
end


function format_dad_iso(PONTOS,SEGMENTOS)
  contorno=calc_ncont(SEGMENTOS);
  ncont=size(contorno,1);
  icrv=0;
  ncurves=sum(contorno[:,2])
  crv=Array{Curve}(undef,ncurves)
  jj=0;
  for j=1:ncont
      pini=contorno[j,1]; # Ponto onde começa o contorno j
      nseg=contorno[j,2];
      for i=1:nseg
          icrv=icrv+1;
          raio=SEGMENTOS[i+pini-1,4] #define valores para o raio
          np1=Int32(SEGMENTOS[i+pini-1,2]); # Define as coordenas x
          np2=Int32(SEGMENTOS[i+pini-1,3]); #Define as coordenasdas y
          p1=[PONTOS[np1,2], PONTOS[np1,3],0]; #Define o primeiro ponto da curva
          p2=[PONTOS[np2,2], PONTOS[np2,3],0]; #Define o segundo ponto da curva
          if(raio==0)
              crv[icrv]=nrbline(p1,p2);
          else
              xc,yc=calcula_centro(p1[1],p1[2],p2[1],p2[2],raio);
              sang,eang = calcula_arco(p1[1],p1[2],p2[1],p2[2],xc,yc,raio);
              centro=[xc,yc,0];
              crv[icrv]=nrbcirc(raio,centro,sang,eang);
          end
      end
  end
  return crv#,contorno
end


function monta_Teq(tipoCDC,valorCDC,x)
# Separa fluxo e temperatura

# ncdc = n�mero de linhas da matriz CDC
# T = vetor que cont�m as temperaturas nos n�s
# q = vetor que cont�m o fluxo nos n�s

ncdc = length(tipoCDC)
T=zeros(ncdc)
q=zeros(ncdc)
for i=1:ncdc # La�o sobre as condi��es de contorno
    if tipoCDC[i] == 1 # Fluxo � conhecido
        T[i] = x[i]; # A temperatura � o valor calculado
        q[i] = valorCDC[i]; # O fluxo � a condi�ao de contorno
    else # A temperatura � conhecida
        T[i] = valorCDC[i]; # A temperatura � a condi�ao de contorno
        q[i] = x[i]; # O fluxo � o valor calculado
    end
end
return T,q
end
function mostra_geo(crvs)
  p=plot(legend=:none,aspect_ratio=:equal)
  for i in crvs
    p=nrbplot(i)
  end
  p
end
