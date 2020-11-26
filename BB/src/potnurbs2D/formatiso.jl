function calcula_arco(x1,y1,x2,y2,xc,yc,raio)
# Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
# horizontal direction 
#  and 
# the tet2 angle between the line from point (x2,y2) to (xc,yc)
# and the horizontal direction
dx1 = x1 - xc; dy1 = y1 - yc
dx2 = x2 - xc; dy2 = y2 - yc
    
if abs(dx1)<1e-15         # <------- EmerEdited
      dx1=0
elseif abs(dy1)<1e-15
        dy1=0
end
# @show dy1,dx1           # <------- EmerEdited

tet1=atan(dy1,dx1);
a=[dx1,dy1,0];
b=[dx2,dy2,0];
c=[0,0,(dx1*dy2-dx2*dy1)];
    
# @show sqrt(sum(c.*c)), sum(a.*b)              # <------- EmerEdited
aa=sqrt(sum(c.*c)); bb=sum(a.*b);              # <------- EmerEdited
if abs(aa)<1e-15                  # <------- EmerEdited
     aa=0
elseif abs(bb)<1e-15                 # <------- EmerEdited
     bb=0
end
    
angle = atan(aa,bb);
    
# angle = atan(norm(cross(a,b)),dot(a,b));  # <------- EmerEdited
# @show tet1,angle           # <------- EmerEdited
if(raio>0)
    tet2=tet1+angle;
    if abs(tet2)<1e-15      # <------- EmerEdited
        tet2=0
    end
else
    tet2=tet1-angle;
    if abs(tet2)<1e-15      # <------- EmerEdited
        tet2=0
    end
end
[tet1,tet2]
end

###############################################################################

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

###############################################################################

function calc_ncont(SEGMENTOS)
# conta o número de contornos e o número de segmentos em cada contorno e 
#     cada linha inicia com o primeiro segmento daquele contorno
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

############################################################################### 

function format_dad_iso(PONTOS,SEGMENTOS) 
#   println("PONTOS = $PONTOS")
#   println("SEGMENTOS = $SEGMENTOS")
  contorno=calc_ncont(SEGMENTOS);
#   println(" contorno = $contorno")
  ncont=size(contorno,1);   # num de contornos
#   @show ncont      # <------------                 emerinserted           emerinserted
#   @show contorno   # <------------                 emerinserted           emerinserted
  
#   println(" ncont = $ncont")
  icrv=0;
#   ncurves=Int64(0);
  ncurves=sum(contorno[:,2])  # num de segmentos ou de curvas
#   println("ncurves = $ncurves")
#   crv=Array{Curve,ncurves} #(ncurves)  # EmerEdited
  crv=Array{Curve}(undef,ncurves)    # EmerEdited
#   jj=0;
  for j=1:ncont  # looping pelos contornos um por um
      pini=contorno[j,1];  # Ponto onde começa o contorno j
      nseg=contorno[j,2];  # num de segmentos do contorno j
      for i=1:nseg   # looping por cada segmento do contorno j
          icrv=icrv+1;
          raio=SEGMENTOS[i+pini-1,4] #define valores para o raio
          np1=Int32(SEGMENTOS[i+pini-1,2]); # Define o primeir o ponto do segmento
          np2=Int32(SEGMENTOS[i+pini-1,3]); #Define o segundo/último ponto do segmento
          p1=[PONTOS[np1,2], PONTOS[np1,3],0]; #Define o primeiro ponto da curva
          p2=[PONTOS[np2,2], PONTOS[np2,3],0]; #Define o segundo/último ponto da curva/segmento
            
          if(raio==0)
              crv[icrv]=nrbline(p1,p2);
          else
              xc,yc=calcula_centro(p1[1],p1[2],p2[1],p2[2],raio);
  #            @show i,p1,p2,xc,yc,raio    #  <-----         emerinserted      emerinserted
              sang,eang = calcula_arco(p1[1],p1[2],p2[1],p2[2],xc,yc,raio);
              centro=[xc,yc,0];

 #             @show sang,eang    #  <-----         emerinserted      emerinserted
              crv[icrv]=nrbcirc(abs(raio),centro,sang,eang); 
          end
      end
  end
  return crv#,contorno
end

###############################################################################

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

###############################################################################

function mostra_geo(crvs)
  p=plot(legend=:none,aspect_ratio=:equal)
  for i in crvs
    p=nrbplot(i)
  end
  p
end
