function gera_vars(ELEM,CCFace,NOS_GEO)

    # Descri��o: Gera a matriz CDC a partir das condi��es de contorno das
    #            faces.
    # Autor:     Gustavo Gontijo
    #

    # Gera a matriz de CDC
    # CDC = [n�mero do elemento, tipo da CDC, valor da CDC no n� 1,...
    #                               valor da CDC no n� 2, valor da CDC no n� 3]
    nelem = size(ELEM,1)
    CDC = zeros(nelem,3)
    NOS=zeros(nelem,4)
    for i=1:nelem
        CDC[i,1] = i
        CDC[i,2] = CCFace[ELEM[i,5],2];
        CDC[i,3] = CCFace[ELEM[i,5],3];
        noselem = ELEM[i,2:4]
        X1=NOS_GEO[noselem[1],2]
        Y1=NOS_GEO[noselem[1],3]
        Z1=NOS_GEO[noselem[1],4]
        X2=NOS_GEO[noselem[2],2]
        Y2=NOS_GEO[noselem[2],3]
        Z2=NOS_GEO[noselem[2],4]
        X3=NOS_GEO[noselem[3],2]
        Y3=NOS_GEO[noselem[3],3]
        Z3=NOS_GEO[noselem[3],4]
        XM=(X1+X2+X3)/3;
        YM=(Y1+Y2+Y3)/3;
        ZM=(Z1+Z2+Z3)/3;
        NOS[i,1]=i
        NOS[i,2]=XM
        NOS[i,3]=YM
        NOS[i,4]=ZM
    end
    return CDC,NOS
end

function calcula_HeGns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsil,w,FR,CW)
    n_pint=length(qsil);
    eta = qsil;
    rho = w;
    Tast = complex(0,0);
    qast = complex(0,0);
    g=0;
    h=0;
    for l=1:n_pint
        for m=1:n_pint
            
            qsi=(1-eta[l])*qsil[m];
            N1,N2,N3=calc_fformatri(qsi,eta[l]);
            x=N1*x1+N2*x2+N3*x3;
            y=N1*y1+N2*y2+N3*y3;
            z=N1*z1+N2*z2+N3*z3;
            
            dNdqsi = [1; 0; -1];
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
            
            Tast,qast=calc_solfund(x,y,z,xd,yd,zd,n,FR,CW);
            
            h=h+qast*(1-eta[l])*rho[l]*w[m]*J;
            g=g+Tast*(1-eta[l])*rho[l]*w[m]*J;
            
        end
    end
    return g,h
end

function calc_fforma_quad(qsi, eta)

    #--------------------------------------------------------------------------
    # Dados de entrada:
    # qsi, eta - pontos de Gauss onde as funções de forma são calculadas.
    # Dados de saída:
    # [N] - Funções de forma para um elemento quadrilateral linear calculadas
    # em (qsi, eta).
    #--------------------------------------------------------------------------

    N = (1. /4.).*[(1.  - qsi).*(1. - eta);
                   (1. + qsi).*(1. - eta);
                   (1. + qsi).*(1. + eta);
                   (1. - qsi).*(1. + eta)]
    N1=N[1]
    N2=N[2]
    N3=N[3]
    N4=N[4]

    return N1, N2, N3, N4

end
# calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,xi,w,k)
function calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,FR,CW)


    # Integra��o singular das matrizes H e G. No caso da matriz G, o elemento
    # triangular � decomposto em tr�s quadril�teros degenerados na forma
    # de tri�ngulos. Os dois primeiros n�s destes quadril�teros s�o
    # coincidentes (formam um s� v�rtice do tri�ngulo) e correspondem ao ponto
    # onde existe a singularidade, ou seja, ao ponto fonte que se localiza no
    # centr�ide do elemento. Isto faz com que haja uma concentra��o de pontos
    # de integra��o junto � singularidade, al�m do jacobiano ser igual a zero
    # na singularidade. No caso da matriz H, a integra��o � anal�tica e sempre
    # ser� igual a -1/2.

    npg=length(qsi); # N�mero de pontos de integra��o
    g = 0.0; # inicializa��o da matriz G

    for kk=1:3
        x1t=xd; # coordenada x do primeiro n� do quadrilatero desgenerado
        y1t=yd; # coordenada y do primeiro n� do quadrilatero desgenerado
        z1t=zd; # coordenada z do primeiro n� do quadrilatero desgenerado
        x2t=xd; # coordenada x do segundo n� do quadrilatero desgenerado
        y2t=yd; # coordenada y do segundo n� do quadrilatero desgenerado
        z2t=zd; # coordenada z do segundo n� do quadrilatero desgenerado
        
        if(kk==1) # Terceiro e quarto n�s do primeiro quadril�tero desgenerado
            x3t=x1; # coordenada x do terceiro n� do quadrilatero desgenerado
            y3t=y1; # coordenada y do terceiro n� do quadrilatero desgenerado
            z3t=z1; # coordenada z do terceiro n� do quadrilatero desgenerado
            x4t=x2; # coordenada x do quarto n� do quadrilatero desgenerado
            y4t=y2; # coordenada y do quarto n� do quadrilatero desgenerado
            z4t=z2; # coordenada z do quarto n� do quadrilatero desgenerado
            
        elseif(kk==2) # Terceiro e quarto n�s do segundo quadril�tero
            # desgenerado
            x3t=x2; # coordenada x do terceiro n� do quadrilatero desgenerado
            y3t=y2; # coordenada y do terceiro n� do quadrilatero desgenerado
            z3t=z2; # coordenada z do terceiro n� do quadrilatero desgenerado
            x4t=x3; # coordenada x do quarto n� do quadrilatero desgenerado
            y4t=y3; # coordenada y do quarto n� do quadrilatero desgenerado
            z4t=z3; # coordenada z do quarto n� do quadrilatero desgenerado
            
        elseif(kk==3) # Terceiro e quarto n�s do terceiro quadril�tero
            # desgenerado
            x3t=x3; # coordenada x do terceiro n� do quadrilatero desgenerado
            y3t=y3; # coordenada y do terceiro n� do quadrilatero desgenerado
            z3t=z3; # coordenada z do terceiro n� do quadrilatero desgenerado
            x4t=x1; # coordenada x do quarto n� do quadrilatero desgenerado
            y4t=y1; # coordenada y do quarto n� do quadrilatero desgenerado
            z4t=z1; # coordenada z do quarto n� do quadrilatero desgenerado
            
        end
        
        for ii = 1: npg # la�o sobre a primeira vari�vel de integra��o
            for jj = 1: npg # la�o sobre a segunda vari�vel de integra��o
                N1,N2,N3,N4 = calc_fforma_quad(qsi[ii],qsi[jj]); # Fun��es de
                #  forma
                xc = x1t*N1+x2t*N2+x3t*N3+x4t*N4; # coordenada x do
                # ponto campo
                yc = y1t*N1+y2t*N2+y3t*N3+y4t*N4; # coordenada y do
                # ponto campo
                zc = z1t*N1+z2t*N2+z3t*N3+z4t*N4; # coordenada z do
                # ponto campo
                J = calc_jacobiano_quad(x1t,y1t,z1t,x2t,y2t,z2t,x3t,y3t,z3t,x4t,y4t,z4t,qsi[ii],qsi[jj]); # jacobiano(varia ao longo
                #  do elemento desgenerado)
                Tast,qast = calc_solfund(xc, yc, zc,xd,yd,zd,[0 0 0], FR,CW) # Sol.
                # fudamental de temperatura
                g = g + w[jj]*w[ii]*J*Tast# integra��o num�rica da matriz G
                
            end
        end
    end

    h=-1/2 # Integra��o anal�tica da matriz H

    return g,h

end

function calc_dfforma_quad(qsi, eta)

    dNdqsi = (1. /4.)*[-(1. - eta);
                       (1. - eta);
                       (1. + eta);
                       -(1. + eta)];

    dNdeta = (1. /4.)*[-(1. - qsi);
                       -(1. + qsi);
                       (1. + qsi);
                       (1. - qsi)];
    return dNdqsi, dNdeta
end

function calc_jacobiano_quad(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,qsi,eta)
    dNdqsi, dNdeta = calc_dfforma_quad(qsi,eta); # Calcula a derivada das funções
    # de forma
    dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3]+x4*dNdqsi[4]
    dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3]+y4*dNdqsi[4]
    dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3]+z4*dNdqsi[4]

    dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3]+x4*dNdeta[4]
    dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3]+y4*dNdeta[4]
    dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3]+z4*dNdeta[4]

    g1 = dydqsi*dzdeta - dzdqsi*dydeta;
    g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
    g3 = dxdqsi*dydeta - dydqsi*dxdeta;
    J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);
    return J
end


function monta_Teq(CDC,x)
    # Separa fluxo e temperatura

    # ncdc = n�mero de linhas da matriz CDC
    # T = vetor que cont�m as temperaturas nos n�s
    # q = vetor que cont�m o fluxo nos n�s

    ncdc = size(CDC,1)
    T=complex(zeros(ncdc))
    q=complex(zeros(ncdc))
    for i=1:ncdc # La�o sobre as condi��es de contorno
        tipoCDC=CDC[i,2]; # Tipo da condi��o de contorno
        valorCDC=complex(CDC[i,3]) # Valor da condi��o de contorno
        valorcalculado=x[i] # Valor que antes era desconhecido
        if tipoCDC == 1 # Fluxo � conhecido
            T[i] = valorcalculado # A temperatura � o valor calculado
            q[i] = valorCDC # O fluxo � a condi�ao de contorno
        else # A temperatura � conhecida
            T[i] = valorCDC # A temperatura � a condi�ao de contorno
            q[i] = valorcalculado # O fluxo � o valor calculado
        end
    end
    return T,q
end

function calc_solfund(x,y,z,xd,yd,zd,n,FR,CW)
    rx=x-xd;
    ry=y-yd;
    rz=z-zd;
    r =sqrt(rx^2+ry^2+rz^2);
    ZW=complex(0.,-FR*r/CW);
    U=exp(ZW)/complex(4*pi*r,0.);
    drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;
    Q=-(ZW-complex(1.,0.))*U*complex(drdn/r,0.);
    return (U),(Q)
end

function calc_fforma(qsi)
    N1=1. /2. .*(1. .-qsi);
    N2=1. /2. .*(1. .+qsi);
    return N1,N2
end


function calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,qsi,eta)
    dNdqsi, dNdeta = calc_dfforma(qsi,eta);
    dxdqsi = x1 .*dNdqsi[1].+ x2.*dNdqsi[2].+x3.*dNdqsi[3]
    dydqsi = y1.*dNdqsi[1].+y2.*dNdqsi[2].+y3.*dNdqsi[3]
    dzdqsi = z1.*dNdqsi[1].+z2.*dNdqsi[2].+z3.*dNdqsi[3]

    dxdeta = x1.*dNdeta[1].+x2.*dNdeta[2].+x3.*dNdeta[3]
    dydeta = y1.*dNdeta[1].+y2.*dNdeta[2].+y3.*dNdeta[3]
    dzdeta = z1.*dNdeta[1].+z2.*dNdeta[2].+z3.*dNdeta[3]

    g1 = dydqsi.*dzdeta .- dzdqsi.*dydeta;
    g2 = dzdqsi.*dxdeta .- dxdqsi.*dzdeta;
    g3 = dxdqsi.*dydeta .- dydqsi.*dxdeta;
    J = sqrt(g1.^2.0 + g2.^2.0 + g3.^2.0);

    return J
end

function calc_fformatri(qsi,eta)
    N1 = qsi;
    N2 = eta;
    N3 = 1 .-qsi .-eta;
    return N1,N2,N3
end
function calc_dfforma(qsi, eta)
    dNdqsi = [1 0 -1]
    dNdeta = [0 1 -1]
    return dNdqsi, dNdeta
end

function calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3)
    v1 = [x3,y3,z3] .- [x2,y2,z2];
    v2 = [x1,y1,z1] .- [x2,y2,z2];
    n = cross(v1, v2);
    if sum(abs.(n)) > 0.000001
        n = n./norm(n);
    end
    return n
end

function Gauss_Legendre(x1,x2,n)
    x = zeros(n);
    w = zeros(n);
    pp=0;
    eps=3e-14;
    m=round(Int,(n+1)/2);
    xm=.5*(x2+x1);
    xl=.5*(x2-x1);
    for i=1:m
        z=cos(pi*(i-.25)/(n+.5));
        while 1==1
            p1=1.;
            p2=0.;
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

function mostra_resultados(XYZ,tri,T)
    nelem=size(tri,1)
    x=zeros(3)
    y=zeros(3)
    z=zeros(3)
    zc=zeros(nelem,1)
    pc=[zeros(3,3)];
    triang=zeros(3,3)
    for elem =1:nelem
        no1=tri[elem,2]
        no2=tri[elem,3]
        no3=tri[elem,4]
        x[1]=XYZ[no1,2]
        y[1]=XYZ[no1,3]
        z[1]=XYZ[no1,4]
        x[2]=XYZ[no2,2]
        y[2]=XYZ[no2,3]
        z[2]=XYZ[no2,4]
        x[3]=XYZ[no3,2]
        y[3]=XYZ[no3,3]
        z[3]=XYZ[no3,4]
        triang=[[x[1] y[1] z[1]
                 x[2] y[2] z[2]
                 x[3] y[3] z[3]]]
        append!(pc,triang)
    end
    fig = plt.figure()
    ax = mp.Axes3D(fig)
    q = ar.Poly3DCollection(pc[2:end], linewidths=1,edgecolors="k")
    ax.add_collection3d(q)
    m = cm.ScalarMappable(cmap=cm.jet)
    b=m.to_rgba(T[1:nelem])
    q.set_facecolor(b[:,1:3])
    m.set_array([minimum(T),maximum(T)])
    m.set_clim(vmin=minimum(T),vmax=maximum(T))
    plt.colorbar(m, orientation="vertical",shrink=0.9)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_xlim(minimum(XYZ[:,2]),maximum(XYZ[:,2]))
    ax.set_ylim(minimum(XYZ[:,3]),maximum(XYZ[:,3]))
    ax.set_zlim(minimum(XYZ[:,4]),maximum(XYZ[:,4]))
    ax.set_aspect("equal")
    ax.view_init(elev=18., azim=43.)
    return
end

function calc_T_pint(PONTOS_int,NOS_GEO,ELEM,T,q,FR,CW,qsi,w,inc)
    # Calcula a temperatura nos pontos internos
    n_pint=length(PONTOS_int[:,1]); # Numero de pontos internos
    n_elem=length(T); # Numero de elementos
    G_int = complex(zeros(n_pint,n_elem))
    H_int = complex(zeros(n_pint,n_elem))
    phi_inc = complex(zeros(n_pint,1))
    g = complex(zeros(n_pint,1))
    for i=1:n_pint # La�o sobre os pontos internos
        x_fonte=PONTOS_int[i,2]; # Coordenada x do ponto fonte
        y_fonte=PONTOS_int[i,3]; # Coordenada y do ponto fonte
        z_fonte=PONTOS_int[i,4]; # Coordenada z do ponto fonte
        for j=1:n_elem  #La�o sobre os elementos
            no1=ELEM[j,2]; # Primeiro n� geom�trico do elemento
            no2=ELEM[j,3]; # Segundo n� geom�trico do elemento
            no3=ELEM[j,4]; # Terceiro n� geom�trico do elemento

            x1=NOS_GEO[no1,2]; # Coordenada x do n� geom�trico 1
            y1=NOS_GEO[no1,3]; # Coordenada y do n� geom�trico 1
            z1=NOS_GEO[no1,4]; # Coordenada z do n� geom�trico 1

            x2=NOS_GEO[no2,2]; # Coordenada x do n� geom�trico 2
            y2=NOS_GEO[no2,3]; # Coordenada y do n� geom�trico 2
            z2=NOS_GEO[no2,4]; # Coordenada x do n� geom�trico 2

            x3=NOS_GEO[no3,2]; # Coordenada x do n� geom�trico 3
            y3=NOS_GEO[no3,3]; # Coordenada y do n� geom�trico 3
            z3=NOS_GEO[no3,4]; # Coordenada z do n� geom�trico 3

            n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio
            G_int[i,j],H_int[i,j]=calcula_HeGns(x1,y1,z1,x2,y2,z2,x3,y3,z3,x_fonte,y_fonte,z_fonte,n,qsi,w,FR,CW); # Chama a functio para calculo de H e G
            # quando o ponto fonte nao pertence ao elemento
        end
        #    if(inc[1,1]!=0)
        #    	#g[i,1]=calc_q(x_fonte,y_fonte,fc,FR,CW,GE);
        #			g[i,1]=0;
        #		else
        #    	g[i,1]=0;
        #    end
        #		if inc[1,1] != 0
	#Vamos incluir um termo de onda incidente
        #				phi_inc[i,1] = calc_inc(x_fonte,y_fonte,z_fonte,FR,CW,inc[1,:]);
        #		end
    end
    T_pint = - (H_int*T - G_int*q - phi_inc)
    #T_pint=-(H_int*T'-G_int*q'-g'); # Vetor que contem a temperatura nos
    #      pontos internos
    return T_pint
end

