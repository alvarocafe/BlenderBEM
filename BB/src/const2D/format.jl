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
    ELEM = zeros(Int,tamanho,3)
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
                xc, yc = calcula_centro(x1l,y1l,x2l,y2l,r);
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
    valorCDC=CCSeg[l,3];
    # valorCDC=complex(CCSeg[l,3],CCSeg[l,4]);
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
return NOS_GEO,NOS,ELEM,CDC,normal
end
function calcula_arco(x1,y1,x2,y2,xc,yc)
    # Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
    # horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
    # and the horizontal direction


    dx1 = x1 - xc; dy1 = y1 - yc;
    dx2 = x2 - xc; dy2 = y2 - yc;

    # Computation of tet1
    if dy1 == 0				# The point 1 and the center have the same y coordinate
        if x1 > xc
            tet1 = 0;
        else  # (x1 < xc)
            tet1 = pi;
        end;
    elseif dx1 == 0				# The point 1 and the center have the same x coordinate
        if y1 > yc
            tet1 = pi/2;
        else  # (y1 < yc)
            tet1 = -pi/2;
        end;
    else  # (dx1~=0 e dy1~=0)
        tet1 = atan(dy1/dx1);
        if dx1<0 && tet1<0
            tet1 = pi + tet1;
        elseif dx1 < 0 && tet1>0
            tet1 = -pi + tet1;
        end;
    end;

    # Computation of tet2
    if dy2 == 0				# The point 2 and the center have the same y coordinate
        if x2 > xc
            tet2 = 0;
        else  # (x2 < xc)
            tet2 = pi;
        end;
    elseif dx2 == 0				# The point 2 and the center have the same x coordinate
        if y2 > yc
            tet2 = pi/2;
        else  # (y2 < yc)
            tet2 = -pi/2;
        end;
    else  # (dx2~=0 e dy2~=0)
        tet2 = atan(dy2/dx2);
        if dx2<0 && tet2<0
            tet2 = pi + tet2;
        elseif dx2 < 0 && tet2>0
            tet2 = -pi + tet2;
        end;
    end;
    return tet1, tet2
end
function calcula_centro(x1,y1,x2,y2,raio)

    # Compute the center of an arc given two points and the radius

    xm=(x1+x2)/2;
    ym=(y1+y2)/2;
    b=sqrt((x2-x1)^2+(y2-y1)^2);
    t1=(x2-x1)/b;
    t2=(y2-y1)/b;
    n1=t2;
    n2=-t1;
    h=sqrt(abs(raio^2-(b/2)^2));
    if(raio>0)
        if(n1==0)
            xc=xm;
            yc=ym-n2/abs(n2)*h;
        else
            xc=-n1/abs(n1)*sqrt(h^2*n1^2/(n1^2+n2^2))+xm;
            yc=n2/n1*(xc-xm)+ym;
        end;
    else
        if(n1==0)
            xc=xm;
            yc=ym+n2/abs(n2)*h;
        else
            xc=n1/abs(n1)*sqrt(h^2*n1^2/(n1^2+n2^2))+xm;
            yc=n2/n1*(xc-xm)+ym;
        end;
    end;
    return xc, yc
end
function calc_fforma(qsi)
    # Evaluates the shape functions for continuous linear elements  
    N1=1/2*(1 .-qsi); # N1 - first shape function for the continuous linear element
    N2=1/2*(1 .+qsi); # N2 - second shape function for the continuous linear element
    return N1,N2
end
function testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio)
    # Rotina para verifica��o da condi��o de um ponto com rela��o
    # � geometria (interno ou externo).
    # O algoritmo utilizado � baseado no n�mero de intersec��es entre
    # o contorno e uma linha vertical que parte do ponto interno.
    # N. �mpar de intersec��es - ponto interno
    # N. par de intersec��es   - ponto externo
    #
    #   Autor: Frederico Lourenao
    #   Data de cria��o setembro de 1999
    #   Revis�o 0.0


    pc = 0;					# Contador de pontos do contorno
    npc = length(xl1);	# N. de pontos do contorno
    sai = false

    while (sai==false && pc < npc)
        pc = pc + 1
        if xpi == xl1[pc]
            xpi = xpi + lx*1e-2
            sai = true
        end
    end

    interv = 0;		# N. de intersec��es entre o contorno e a linha
    # vertical abaixo do ponto interno.
    l=1
    while l<= npc; # for over the lines that defines the boundary
        x1 = xl1[l]; y1 = yl1[l];
        x2 = xl2[l]; y2 = yl2[l]
        if(raio[l]==0) # The segment is a straight line
            if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1) # Check only if the x coordinate of the point is
                # between the x coordinate of the intial and final points of the recta
                m = (y2-y1)/(x2-x1);		# Angular coeficient of the recta
                yi = m*(xpi-x1)+y1;			# y coordinate in the intersection
                if ypi  >yi # If this is true; the point is above the line
                    interv = interv + 1;	# Counter the intersection
                end
            end
        else   # The segment is an arc
            xc,yc=calcula_centro(x1,y1,x2,y2,raio[l]); # compute the center of the arc
            if(xpi<xc+abs(raio[l]) && xpi>xc-abs(raio[l])) # check only if the x coordinate of the point is between
                # the values xc-radius and xc+radius
                teta_i,teta_f=calcula_arco(x1,y1,x2,y2,xc,yc); # compute the arc between the line that defines the arc and
                # the horizontal direction (-π<teta<π)
                tetac=zeros(2,1)
                y=zeros(2,1)
                tetac[1]=acos((xpi-xc)/abs(raio[l])); # first intersection of the horizontal line that cross the
                # point with the circunference (angle between the radius that cross the intersection point and
                # the horizontal direction)
                tetac[2]=-tetac[1];  # angle of the second intersection
                y[1]=abs(raio[l])*sin(tetac[1])+yc; # y coordinate of the first intersection point
                y[2]=abs(raio[l])*sin(tetac[2])+yc; # y coordinate of the second intersection point
                for k=1:2 # check if the angles of the two intersection points are between the angle of the initial and the
                    # final angles of the arc. If yes so the vertical line from the the point to the botton direction
                    # intercept the arc (the counter should be increased).
                    if(raio[l]>0) # the center is on the left of the arc (from the initial point to the final point)
                        if(teta_f>teta_i)
                            if(tetac[k]>teta_i && tetac[k]<teta_f)
                                if(y[k]<ypi)
                                    interv=interv+1
                                end
                            end
                        else # teta_f<teta_i
                            if(tetac[k]>teta_i || tetac[k]<teta_f)
                                if(y[k]<ypi)
                                    interv=interv+1
                                end
                            end
                        end
                    else # raio[l] < 0 the center is on the right of the arc (from the initial point to the final point)
                        if(teta_i > teta_f)
                            if(tetac[k]>teta_f && tetac[k]<teta_i)
                                if(y[k]<ypi)
                                    interv=interv+1
                                end
                            end
                        else # teta_i < teta_f
                            if(tetac[k]>teta_f || tetac[k]<teta_i)
                                if(y[k]<ypi)
                                    interv=interv+1
                                end
                            end
                        end
                    end
                end
            end
        end
        l=l+1
    end


    if rem(interv,2) != 0	# Resto da divis�o de interv por 2
        ponto = "interno"
    else
        ponto = "externo"
    end
    return  xpi,ponto
end

function  aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio)
    # Fun??o para testar a proximidade de um ponto interno ao contorno
    # Se o ponto estiver mais proximo que o desejado (d<d_min) o ponto
    # ? reprovado (aceita = nao). Caso contrr?rio (aceita = sim)


    pc = 0;					# Contador de pontos do contorno
    npc = length(xl1);	# N. de pontos do contorno
    aceita = "sim"

    # Verificando a proximidade do ponto interno aos pontos do contorno
    while (aceita=="sim"  && pc < npc)
        pc = pc + 1
        x1 = xl1[pc]
        y1 = yl1[pc]
        d = √((xpi-x1)^2+(ypi-y1)^2)
        if d < d_min
            aceita = "nao"
        end
    end

    # Verificando a proximidade ?s linhas do contorno
    l = 0;	# Contador das linhas do contorno

    while (aceita=="sim" && l < npc)
        l = l + 1
        x1 = xl1[l]; y1 = yl1[l]
        x2 = xl2[l]; y2 = yl2[l]
        if(raio[l]==0) # The segment is a straight line

            if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1)
                m = (y2-y1)/(x2-x1);		# Angular coeficient of the recta
                yi = m*(xpi-x1)+y1;			# y coordinate of the intersection of a recta normal to the current boundary
                # recta that cross the point
                dy = ypi - yi
                d = abs(dy*cos(atan(m))); # distance from the point to the recta
                if d < d_min
                    aceita = "nao"
                end
            end

            if (ypi > y1 && ypi < y2) || (ypi > y2 && ypi < y1)
                if x1 == x2
                    d = abs(xpi-x1)
                else
                    m = (y2-y1)/(x2-x1);			# Angular coeficient of the recta
                    xi = 1/m*(ypi-y1)+x1;	# x coordinate of the intersection of a recta normal to the current boundary
                    # recta that cross the point
                    dx = xpi - xi
                    d = abs(dx*sin(atan(m))); # distance from the point to the recta
                end
                if d < d_min
                    aceita = "nao"
                end
            end
        else    # The segment is an arc
            xc,yc=calcula_centro(x1,y1,x2,y2,raio[l]); # Center of the arc
            teta_i,teta_f=calcula_arco(x1,y1,x2,y2,xc,yc); # angle of the lines that defines the arc with the horizontal direction
            teta_i,teta_p=calcula_arco(x1,y1,xpi,ypi,xc,yc);  # teta_p angle of the line that cross the center point and the
            # internal point with the horizontal direction
            if(raio[l]>0) # The center is in the left side of the arc (from the initial to the end point)
                if(teta_f>teta_i)
                    if(teta_p>teta_i && teta_p<teta_f)
                        d=abs(raio[l]-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
                        if d < d_min
                            aceita = "nao"
                        end
                    end
                else # teta_f<teta_i
                    if(teta_p>teta_i || teta_p<teta_f)
                        d=abs(raio[l]-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
                        if d < d_min
                            aceita = "nao"
                        end
                    end
                end
            else # raio[l] < 0 # The center is in the right side of the arc (from the initial to the end point)
                if(teta_i > teta_f)
                    if(convert(Float64,teta_p)>convert(Float64,teta_f) && convert(Float64,teta_p)<convert(Float64,teta_i))
                        d=abs(abs(raio[l])-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
                        if d < d_min
                            aceita = "nao"
                        end
                    end
                else # teta_i < teta_f
                    if(teta_p>teta_f || teta_p<teta_i)
                        d=abs(abs(raio[l])-√((xpi-xc)^2+(ypi-yc)^2)); # distance from the point to the arc
                        if d < d_min
                            aceita = "nao"
                        end
                    end
                end
            end
        end
    end
    # Retorna a aprova??o ou nao do ponto interno
    return aceita
end

function gera_p_in(NPX,NPY,PONTO,SEGMENTOS)
    #PONTOS_INT
    # Programa para cria��o de pontos internos a uma  geometria
    # gen�rica formada por retas
    #
    #   Autor: Frederico Lourenao
    #   Data de cria��o setembro de 1999
    #   Revis�o 0.0

    # Defini��o da �rea m�xima para cria��o de pontos internos
    xmin = minimum(PONTO[:,2])
    xmax = maximum(PONTO[:,2])
    ymin = minimum(PONTO[:,3])
    ymax = maximum(PONTO[:,3])
    lx = xmax - xmin;		# Largura do ret�ngulo que cont�m a geometria
    ly = ymax - ymin;		# Altura do ret�ngulo que cont�m a geometria
    n_SEGMENTOSs = length(SEGMENTOS[:,1])

    # Defini��o da maior SEGMENTOS do problema
    max_dl = 0
    for lin = 1 : length(SEGMENTOS[:,1])
        p1 = round(UInt64,SEGMENTOS[lin,2])
        p2 = round(UInt64,SEGMENTOS[lin,3])
        xp1 = PONTO[p1,2]
        yp1 = PONTO[p1,3]
        xp2 = PONTO[p2,2]
        yp2 = PONTO[p2,3]
        dl = √((xp1-xp2)^2+(yp1-yp2)^2)
        if dl > max_dl
            max_dl = dl
        end
    end


    d_min = 0.003*max_dl;	# Dist�ncia m�nima dos pontos internos ao contorno
    npx = NPX+1;				# N. de pontos na horizontal
    npy = NPY+1;				# N. de pontos na vertical

    PONTOS_INT =zeros(NPX*NPY,2)
    # Atribui��o dos pontos finais e iniciais das SEGMENTOSs aos
    # vetores xl1; xl2; yl1 e yl2
    xl1=zeros(n_SEGMENTOSs)
    xl2=zeros(n_SEGMENTOSs)
    yl1=zeros(n_SEGMENTOSs)
    yl2=zeros(n_SEGMENTOSs)
    raio=zeros(n_SEGMENTOSs)

    for t = 1 : n_SEGMENTOSs		# Percorre todas as SEGMENTOSs
        xl1[t] = PONTO[round(UInt64,SEGMENTOS[t,2]),2]
        xl2[t] = PONTO[round(UInt64,SEGMENTOS[t,3]),2]
        yl1[t] = PONTO[round(UInt64,SEGMENTOS[t,2]),3]
        yl2[t] = PONTO[round(UInt64,SEGMENTOS[t,3]),3]
        raio[t]= SEGMENTOS[t,4]
    end

    npi = 0;	# Inicializa��o no n. de pontos internos
    for i = 1 : NPY
        # Cria��o do candidato a ponto interno (xpi,ypi)
        ypi = ymin + (ly/npy)*i;	# y dentro do ret�ngulo
        for j = 1 : NPX
            xpi = xmin + (lx/npx)*j;	# x dentro do ret�ngulo

            # In�cio dos testes para valida��o do ponto interno

            # 1. Verificando se o ponto est� dentro da geometria
            xpi,ponto = testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio)
            # 2. Verificando se o ponto est� muito pr�ximo do contorno
            if (ponto=="interno")
                aceita = aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio)
            else
                aceita = "nao"
            end
            # Armazenando os dados do ponto interno
            if (aceita=="sim")	# O ponto est� dentro da geometria e
                npi = npi + 1;		# e respeita a dist�ncia ao contorno
                PONTOS_INT[npi,:] = [xpi ypi]

            end
        end
    end
    PONTOS_INT=PONTOS_INT[1:npi,:]
    return  PONTOS_INT
end
# Nesse ponto est�o calculados os pontos internos
function inpoly(p, node, edge, reltol = 1.0e-12)
  
    #   INPOLY: Point-in-polygon testing.
    # 
    #  Determine whether a series of points lie within the bounds of a polygon
    #  in the 2D plane. General non-convex, multiply-connected polygonal
    #  regions can be handled.
    # 
    #  SHORT SYNTAX:
    # 
    #    in = inpoly(p, node);
    # 
    #    p   : The points to be tested as an Nx2 array [x1 y1; x2 y2; etc].
    #    node: The vertices of the polygon as an Mx2 array [X1 Y1; X2 Y2; etc].
    #          The standard syntax assumes that the vertices are specified in
    #          consecutive order.
    # 
    #    in  : An Nx1 logical array with IN(i) = TRUE if P(i,:) lies within the
    #          region.
    # 
    #  LONG SYNTAX:
    # 
    #   [in, on] = inpoly(p, node, edge, tol);
    # 
    #   edge: An Mx2 array of polygon edges, specified as connections between
    #         the vertices in NODE: [n1 n2; n3 n4; etc]. The vertices in NODE
    #         do not need to be specified in connsecutive order when using the
    #         extended syntax.
    # 
    #   on  : An Nx1 logical array with ON(i) = TRUE if P(i,:) lies on a
    #         polygon edge. (A tolerance is used to deal with numerical
    #         precision, so that points within a distance of
    #         reltol*min(bbox(node)) from a polygon edge are considered "on" the 
    #         edge.
    # 
    #  EXAMPLE:
    # 
    #    polydemo;       #  Will run a few examples
    # 
    #  See also INPOLYGON
  
    #  The algorithm is based on the crossing number test, which counts the
    #  number of times a line that extends from each point past the right-most
    #  region of the polygon intersects with a polygon edge. Points with odd
    #  counts are inside. A simple implementation of this method requires each
    #  wall intersection be checked for each point, resulting in an O(N*M)
    #  operation count.
    # 
    #  This implementation does better in 2 ways:
    # 
    #    1. The test points are sorted by y-value and a binary search is used to
    #       find the first point in the list that has a chance of intersecting
    #       with a given wall. The sorted list is also used to determine when we
    #       have reached the last point in the list that has a chance of
    #       intersection. This means that in general only a small portion of
    #       points are checked for each wall, rather than the whole set.
    # 
    #    2. The intersection test is simplified by first checking against the
    #       bounding box for a given wall segment. Checking against the bbox is
    #       an inexpensive alternative to the full intersection test and allows
    #       us to take a number of shortcuts, minimising the number of times the
    #       full test needs to be done.
    # 
    #    Darren Engwirda: 2005-2009
    #    Email          : d_engwirda@hotmail.com
    #    Last updated   : 28/03/2009 with MATLAB 7.0
    # 
    #  Problems or suggestions? Email me.
    nnode = size(node,1);
  
    # #  PRE-PROCESSING
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    n  = size(p,1);
    nc = size(edge,1);
  
    #  Choose the direction with the biggest range as the "y-coordinate" for the
    #  test. This should ensure that the sorting is done along the best
    #  direction for long and skinny problems wrt either the x or y axes.
    dxy = maximum(p,1)-minimum(p,1);
    if dxy[1]>dxy[2]
        #  Flip co-ords if x range is bigger
        p = p[:,[2,1]];
        node = node[:,[2,1]];
    end
  
    #  Polygon bounding-box
    dxy = maximum(node,1)-minimum(node,1);
    tol = reltol*minimum(dxy);
    if tol==0.0
        tol = reltol;
    end
  
    #  Sort test points by y-value
    i = sortperm(p[:,2]);
    y= p[i,2];
    x = p[i,1];
    # #  MAIN LOOP
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    cn = falses(n);     #  Because we're dealing with mod(cn,2) we don't have
    #  to actually increment the crossing number, we can
    #  just flip a logical at each intersection (faster!)
    on = cn[:];
    for k = 1:nc         #  Loop through edges
        #  Nodes in current edge
        n1 = edge[k,1];
        n2 = edge[k,2];
      
        #  Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
        #            - also get xmin = min(x1,x2), xmax = max(x1,x2)
        y1 = node[n1,2];
        y2 = node[n2,2];
        if y1<y2
            x1 = node[n1,1];
            x2 = node[n2,1];
        else
            yt = y1;
            y1 = y2;
            y2 = yt;
            x1 = node[n2,1];
            x2 = node[n1,1];
        end
        if x1>x2
            xmin = x2;
            xmax = x1;
        else
            xmin = x1;
            xmax = x2;
        end
        #  Binary search to find first point with y<=y1 for current edge
        if y[1]>=y1
            start = 1;
        elseif y[n]<y1
            start = n+1;       
        else
            lower = 1;
            upper = n;
            for j = 1:n
                start = convert(Int32,round(1/2*(lower+upper)));
                if y[start]<y1
                    lower = start+1;
                elseif y[start-1]<y1
                    break;
                else
                    upper = start-1;
                end
            end
        end
        #  Loop through points
        for j = start:n
            #  Check the bounding-box for the edge before doing the intersection
            #  test. Take shortcuts wherever possible!
            Y = y[j];   #  Do the array look-up once & make a temp scalar
            if Y<=y2
                X = x[j];   #  Do the array look-up once & make a temp scalar
                if X>=xmin
                    if X<=xmax
                      
                        #  Check if we're "on" the edge
                        on[j] = on[j] || (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<=tol);
                      
                        #  Do the actual intersection test
                        if (Y<y2) && ((y2-y1)*(X-x1)<(Y-y1)*(x2-x1))
                            cn[j] = ~cn[j];
                        end
                      
                    end
                elseif Y<y2   #  Deal with points exactly at vertices
                    #  Has to cross edge
                    cn[j] = ~cn[j];
                end
            else
                #  Due to the sorting, no points with >y
                #  value need to be checked
                break
            end
        end
    end
#  Re-index to undo the sorting
cn[i] = cn.|on;
on[i] = on;

return cn;

end      #  inpoly()
function mostra_problema(ELEM,NOS_GEO,NOS,tipoCDC,valorCDC,normal,T,q)
    nelem=size(ELEM,1);
    figure();
    ind=collect(1:nelem);
    indKnownq=ind[tipoCDC];
    indKnownT=ind[.!tipoCDC];
    maxT=maximum(abs.(T[indKnownT]))
    maxq=maximum(abs.(q[indKnownq]))
    xmax=maximum(NOS_GEO[:,1])
    ymax=maximum(NOS_GEO[:,2])
    xmin=minimum(NOS_GEO[:,1])
    ymin=minimum(NOS_GEO[:,2])
    deltax=xmax-xmin
    deltay=ymax-ymin
    dmax=sqrt(deltax^2+deltay^2)
    ax = plt.gca() # get current axes
    plt.grid("on")
    plt.quiver(NOS[indKnownT,1],NOS[indKnownT,2],normal[1,indKnownT].*valorCDC[indKnownT],
               normal[2,indKnownT].*valorCDC[indKnownT],color="red",width=0.002,scale=50, headaxislength=0)
    plt.quiver(NOS[indKnownq,1],NOS[indKnownq,2],normal[1,indKnownq].*valorCDC[indKnownq],
               normal[2,indKnownq].*valorCDC[indKnownq],color="blue",width=0.002,scale=50)
    plt.plot(NOS[indKnownT,1],NOS[indKnownT,2],"ro",markersize=4);	# Plot the node of the elements
    plt.plot(NOS[indKnownq,1],NOS[indKnownq,2],"bo",markersize=4);	# Plot the node of the elements
    ELEM2=[ELEM ELEM[:,2]];
    plt.triplot(NOS_GEO[:,1], NOS_GEO[:,2], ELEM2-1, color=(0.0,0.,0.),linewidth=0.4) 	# Plot the countour of the problem
    plt.axis("equal")
    ax[:set_xlim]((xmin-.15*deltax,xmax+.15*deltax));
    ax[:set_ylim]((ymin-.15*deltay,ymax+.15*deltay));
end

function mostra_heatmap(NOS,PONTOS_INT,T,Ti,NOS_GEO,ELEM,dTdx,dTdy)
    nelem=size(ELEM,1);
    figure();
    xmax=maximum(NOS_GEO[:,1])
    ymax=maximum(NOS_GEO[:,2])
    xmin=minimum(NOS_GEO[:,1])
    ymin=minimum(NOS_GEO[:,2])
    deltax=xmax-xmin
    deltay=ymax-ymin
    dmax=sqrt(deltax^2+deltay^2)
  
    XY=[NOS;PONTOS_INT];
    triang = tri.Triangulation(XY[:,1], XY[:,2])
    cor=[T; Ti]
    t=triang[:triangles]+1
    centroid=(XY[t[:,1],:]+XY[t[:,2],:]+XY[t[:,3],:])/3.0;
    i = inpoly(centroid,NOS,ELEM);
    ind=collect(1:size(t,1));
    ind=ind[i];                                   
    #  Take triangles with internal centroids
    t = t[ind,:];
  
    nelem=size(ELEM,1);
    plt.plot(NOS_GEO[:,1],NOS_GEO[:,2],"kx",markersize=4,linewidth=1);	# Plot the node of the elements
    plt.grid("on")
  
    plt.plot(NOS[:,1],NOS[:,2],"ko",markersize=4);	# Plot the node of the elements
    ELEM2=[ELEM ELEM[:,2]];
  
    ncont=50;
    # plt.triplot(XY[:,1], XY[:,2], t-1, color=(0.0,0.,0.),linewidth=0.4)
    plt.triplot(NOS_GEO[:,1], NOS_GEO[:,2], ELEM2-1, color=(0.0,0.,0.),linewidth=0.4)
    plt.tricontourf(XY[:,1], XY[:,2], t-1, cor, ncont)
    plt.colorbar()
  
    plt.quiver(PONTOS_INT[:,1],PONTOS_INT[:,2],dTdx,dTdy,color="red",width=0.002,scale=50)
    plt.axis("equal")
    ax = plt.gca() # get current axes
    ax[:set_xlim]((xmin-.15*deltax,xmax+.15*deltax));
    ax[:set_ylim]((ymin-.15*deltay,ymax+.15*deltay));
  
end
function telles(gamm,eet)

    eest = eet^2 - 1;
    term1 = eet*eest + abs(eest);
    if term1 < 0
        term1 = (-term1)^(1/3);
        term1 = -term1;
    else
        term1 = term1^(1/3);
    end

    term2 = eet*eest - abs(eest);
    if term2 < 0
        term2 = (-term2)^(1/3);
        term2 = -term2;
    else
        term2 = term2^(1/3);
    end
    GAMM = term1 + term2 + eet;


    Q = 1 + 3*GAMM^2;
    A = 1/Q;
    B = -3*GAMM/Q;
    C = 3*GAMM^2/Q;
    D = -B;

    eta = A*gamm.^3 .+ B*gamm.^2 .+ C*gamm .+ D;
    Jt = 3*A*gamm.^2 .+ 2*B*gamm .+ C;
    return  eta,Jt
end
function trocacol(G,H,CDC)
    # Applies the boundary conditions by exchanging the columns of matrices H and G
    ncdc = length(CDC); # Number of boundary conditions
    A=H*1.0;	# Make a copy of H
    B=G*1.0;	# Make a copy of G
    for i=1:ncdc # Loop over the boundary conditions
        tipoCDC = CDC[i]; # Type of boundary condition
        if tipoCDC == 0 # The velocity potential is known
            colunaA=-A[:,i]; # 
            A[:,i]=-B[:,i]; # Matrix A receives the column from matrix G
            B[:,i]=colunaA; # Matrix B receives the column from matrix H
        end
    end
    #  b=B*valoresconhecidos; # Vector b
    return A, B
end
function Gauss_Legendre(x1,x2,n)
  eps=3e-14;
  m::Int64 = round((n+1)/2);
  x = zeros(n)
  w = zeros(n)
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
