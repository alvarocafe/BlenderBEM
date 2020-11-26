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

function format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg)
  #Programa para formatar os dados de entrada
  #Autor: Álvaro Campos Ferreira, Éder Lima de Albuquerque

  #Vamos declarar as variáveis que estão com type instabilities
  r1::Float64 = 0;
  r2::Float64 = 0;
  tet::Float64 = 0;
  tet2::Float64 = 0;
  sig::Int64 = 1;
  divtet::Float64 = 0;
  x_i::Float64 = 0;
  y_i::Float64 = 0;
  x_m::Float64 = 0;
  y_m::Float64 = 0;
  x_f::Float64 = 0;
  y_f::Float64 = 0;

  #Definimos o tamanho das matrizes
  tamanho::Int64 = 0   #Iniciamos o tamanho

  for i = 1 : size(MALHA)[1]
    tamanho = tamanho + MALHA[i,2]
  end
  NOS_GEO = zeros(tamanho,3)
  ELEM = zeros(tamanho,3)
  NOS = zeros(tamanho,3)
  normal = zeros(tamanho,tamanho)
  cont_nos = 0;
  cont_el = 0;
  num_lin = size(SEGMENTOS)[1]
  p_ini::Int64 = SEGMENTOS[1,2]
  #Definição da maior dimensão do problema
  max_dl::Float64 = 0
  for lin = 1 : num_lin
    p1::Int64 = SEGMENTOS[lin,2]
    p2::Int64 = SEGMENTOS[lin,3]
    xp1 = PONTOS[p1,2]
    yp1 = PONTOS[p1,3]
    xp2 = PONTOS[p2,2]
    yp2 = PONTOS[p2,3]
    dl = sqrt((xp1 - xp2)^2 + (yp1 - yp2)^2)
    if dl > max_dl
      max_dl = dl
    end
  end

  no_ini = 1
  t = 1
  pp2::Int64 = 2
  no1_prox = 0
  while(t<num_lin)  	# While over all lines
    while(pp2!=p_ini)
      num_el_lin = MALHA[t,2];	# Number of the elements in the line t
      # Coordinates of the initial and final PONTOS of each line
      # (x1l,y1l,x2l,y2l)
      pp1::Int64  = SEGMENTOS[t,2];
      pp2 = SEGMENTOS[t,3];
      x1l = PONTOS[pp1,2];
      y1l = PONTOS[pp1,3];
      x2l = PONTOS[pp2,2];
      y2l = PONTOS[pp2,3];
      # 1. Generation of the matrices NOS, NOS_GEO, NOS_DRM, ELEM e ELEM_GEO
      if(SEGMENTOS[t,4]==0) # The segment is a straight line
        # Increment in x and y direction
        delta_x = x2l - x1l;
        delta_y = y2l - y1l;
      else #The segment is an arc
        # Compute the center of the arc and its coordinates
        r = SEGMENTOS[t,4];
        xc,yc = calcula_centro(x1l,y1l,x2l,y2l,r);
        # Distance between p1 and c (r1) and between p2 and c (r2)
        r1 = sqrt((x1l-xc)^2+(y1l-yc)^2);
        r2 = sqrt((x2l-xc)^2+(y2l-yc)^2);
        if abs(r1-r2)<.00001*max_dl
          # Compute the angle between the lines from point c to p1 (tet1) and c to p2 (tet2)
          tet1, tet2 = calcula_arco(x1l,y1l,x2l,y2l,xc,yc);
          if tet2 < tet1
            tet2 = tet2 + 2*pi;
          end;

          # Angle of the sector defined by the arc
          if SEGMENTOS[t,4] > 0
            tet = abs(tet2-tet1);
            sig = 1;
          else
            tet = 2*pi-abs(tet2-tet1);
            sig = -1;
          end;

          # Angle between two nodes of the line
          divtet = tet/(2*num_el_lin);
        else
          println("Error in the data input file: Wrong central point");
        end;
      end
      # Generation of elements and nodes
      for i = 1 : num_el_lin

        if(SEGMENTOS[t,4]==0) # The segment is a straight line
          x_i = x1l + delta_x/num_el_lin*(i-1);			# initial x coordinate of the element
          y_i = y1l + delta_y/num_el_lin*(i-1);			# initial y coordinate of the element
          x_m = x1l + delta_x/num_el_lin*(i-.5);	# midpoint x coordinate of the element
          y_m = y1l + delta_y/num_el_lin*(i-.5);	# midpoint y coordinate of the element
          lx=x_m-x_i;                              # distance in x direction between geometric nodes 1 and 2
          ly=y_m-y_i;                              # distance in y direction between geometirc nodes 1 and 2

        else  # The segment is an arc
          # Compute the node coordinates
          x_i = xc+r1*cos(tet1+2*(i-1)*sig*divtet);
          y_i = yc+r1*sin(tet1+2*(i-1)*sig*divtet);
          x_m = xc+r1*cos(tet1+(2*i-1)*sig*divtet);
          y_m = yc+r1*sin(tet1+(2*i-1)*sig*divtet);
        end;
        cont_el = cont_el + 1;
        # Set the coordinate of the physical and DRM nodes
        if(no1_prox==0) # the first node needs to be created
          cont_nos = cont_nos + 1;
          NOS_GEO[cont_nos,:]=[cont_nos,x_i,y_i];
          no1=cont_nos;
        else
          no1=no1_prox;
        end;
        if(pp2!=p_ini || i<num_el_lin)
          cont_nos = cont_nos + 1;
          if(SEGMENTOS[t,4]==0) # Straight line
            x_f=x_m+lx;
            y_f=y_m+ly;
          else                 # arc
            x_f = xc+r1*cos(tet1+(2*i)*sig*divtet);
            y_f = yc+r1*sin(tet1+(2*i)*sig*divtet);
          end;
          NOS_GEO[cont_nos,:]=[cont_nos,x_f,y_f];
          no3=cont_nos;
          no1_prox=no3;
        else
          if(no_ini==0)
            cont_nos = cont_nos + 1;
            if(SEGMENTOS[t,4]==0) # Straight line
              x_f=x_m+lx;
              y_f=y_m+ly;
            else                 # Arc
              x_f = xc+r1*cos(tet1+(2*i)*sig*divtet);
              y_f = yc+r1*sin(tet1+(2*i)*sig*divtet);
            end;
            NOS_GEO[cont_nos,:]=[cont_nos,x_f,y_f];
            no3=cont_nos;
            no1_prox=0;
          else
            no3=no_ini;
            no1_prox=0;
          end;
        end;
        ELEM[cont_el,:]=[cont_el,no1,no3];
      end;# end of for i = 1 : num_el_lin
      if pp2 == p_ini
        if t < num_lin
          p_ini = SEGMENTOS[t+1,2];
          if(SEGMENTOS[t+1,3]==2)
            no_ini = 0;
          else
            no_ini = cont_nos+1;
          end
        end;
      end;
      t=t+1;
    end;                                  #end of while p2
  end

  # Gera��o da matriz CDC (Condi��es de Contorno)
  # CDC = [n. do elemento, tipo de cdc, valor da cdc]
  # Tipos de cdc: 0 : temperatura conhecida
  #               1 : fluxo conhecido
  n_elem=size(ELEM,1);
  cont_el2 = 0;
  CDC=zeros(n_elem,3);
  for l = 1 : length(SEGMENTOS[:,1])
    n_el_lin = MALHA[l,2];
    el_ini::Int64 = cont_el2 + 1;
    el_fin::Int64 = cont_el2 + n_el_lin;
    valorCDC=CCSeg[l,3]
    tipoCDC=CCSeg[l,2];
    for el = el_ini : el_fin
      CDC[el,:] = [el,tipoCDC,valorCDC];
    end;
    cont_el2 = el_fin;
  end;

  nnos=size(ELEM,1); # N�mero de n�s
  NOS=zeros(nnos,3);
  for i=1:nnos # La�o sobre os pontos fontes
    pontoi::Int64=ELEM[i,2]; # Ponto final do elemento
    pontof::Int64=ELEM[i,3]; # Ponto inicial do elemento
    xi=NOS_GEO[pontoi,2]; # Coordenada x do ponto inicial do elemento
    xf=NOS_GEO[pontof,2]; # Coordenada x do ponto final do elemento
    yi=NOS_GEO[pontoi,3]; # Coordenada y do ponto inicial do elemento
    yf=NOS_GEO[pontof,3];  # Coordenada y do ponto final do elemento
    xd=(xi+xf)/2; # Coordenada x do ponto fonte
    yd=(yi+yf)/2; # Coordenada y do ponto fonte
    h1=xf-xi;
    h2=yf-yi;
    el = sqrt(h1^2 + h2^2);
    normal[1,i] = h2/el;
    normal[2,i] = -h1/el;
    NOS[i,1]=i;
    NOS[i,2]=xd;
    NOS[i,3]=yd;
  end
  return NOS_GEO,NOS,ELEM,CDC
end


