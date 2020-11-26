function cal_Aeb(b1,b2,collocCoord,nnos,crv,dcrv,E,kmat,CDC)
   # Solves the Helmholtz equation
    n = length(crv);	# Number of curves
    ncollocpoints=size(collocCoord,1)
    H=complex(zeros(ncollocpoints,ncollocpoints));
    G=complex(zeros(ncollocpoints,ncollocpoints));
    npgauss=12;
    qsi,w=Gauss_Legendre(-1,1,npgauss)
    # qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1 : n #Laço sobre as curvas (elementos)
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
    H=H+E/2;
    # Aplica as condições de contorno trocando as colunas das matrizes H e G
    ncdc = size(CDC,1); # número de linhas da matriz CDC
    A=copy(H);
    B=copy(G);
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

function cal_Aebpot(b1,b2,collocCoord,nnos,crv,dcrv,E,kmat,CDC)
   # Solves the Helmholtz equation
    n = length(crv);	# Number of curves
    ncollocpoints=size(collocCoord,1)
    H=complex(zeros(ncollocpoints,ncollocpoints));
    G=complex(zeros(ncollocpoints,ncollocpoints));
    npgauss=12;
    qsi,w=Gauss_Legendre(-1,1,npgauss)
    # qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1 : n #Laço sobre as curvas (elementos)
        uk=unique(crv[i].knots)
        ne=length(uk)-1;
        for j=1:ne
            range=uk[j+1]-uk[j]
            mid=(uk[j+1]+uk[j])/2
            for k=1:ncollocpoints
                xfonte=collocCoord[k,1]
                yfonte=collocCoord[k,2]
                # Integrais de dom�nio
                g,h,id=integra_elem_POT(xfonte,yfonte,crv[i],dcrv[i],range,mid,qsi,w,kmat); # Integra��o sobre o
                #  element (I = int F n.r/r dGama)
                H[k,id+nnos[i]]=H[k,id+nnos[i]]+h;
                G[k,id+nnos[i]]=G[k,id+nnos[i]]+g;
            end
        end
    end
    H=H-E/2;
    # Aplica as condições de contorno trocando as colunas das matrizes H e G
    ncdc = size(CDC,1); # número de linhas da matriz CDC
    A=copy(H);
    B=copy(G);
    for i=1 : ncdc # Laço sobre as condições de contorno
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


function calcula_iso(collocCoord,nnos,crv,dcrv,E,kmat)
    # Solves the Helmholtz equation
    n = length(crv);	# Number of curves
    ncollocpoints=size(collocCoord,1)

    H=complex(zeros(ncollocpoints,ncollocpoints));
    G=complex(zeros(ncollocpoints,ncollocpoints));
    npgauss=12;
    qsi,w=Gauss_Legendre(-1,1,npgauss)
    # qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1 : n #Laço sobre as curvas (elementos)
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
    H=H+E/2;
    return H,G
end

function  integra_elem(xfonte,yfonte,crv,dcrv,range,mid,qsi,w,k)
    # Helmholtz equation
    # Integra��o sobre os elementos (integral I)
    # N�mero de pontos de Gauss usados na integra��o de I
    npg = length(qsi);
    dudqsi=range/2
    g = complex(zeros(crv.order))
    h = complex(zeros(crv.order))
    id=zeros(Int,crv.order,1)
    # dudqsi=(range[2]-range[1]) / 2;
    # eet = (range[1]+range[2]-2*cpt)/(range[1]-range[2]);
    # eta,Jt=telles(qsi,eet);
    for i = 1 : npg # Percorre os pontos de integra��o
        # qsi_param=convertToParamSpace(eta[i], range);
        p,dp=nrbdeval(crv, dcrv, qsi[i]/2*range+mid)
        B, id = nrbbasisfun(crv,qsi[i]/2*range+mid)
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
        
        ZR=real(k*r);
        Z=complex(0.,ZR);
        F0C=besselk(0,Z);
        F1C=besselk(1,Z);
        qast=-(Z/r*nr*F1C)/(2*pi); #Solução Fundamental da pressão acústica
        Tast=F0C/(2*pi);    #Solução Fundamental do fluxo de pressão acústica
        #qast = complex(1,0)
        #Tast = complex(1,0)
        h = h + B*qast*dgamadu *dudqsi* w[i]
        g = g + B*Tast*dgamadu *dudqsi* w[i]
    end
    
    return g,h,id
end

function calcula_iso_POT(collocCoord,nnos,crv,dcrv,E,kmat)
    # Solves the Laplace equation


    n = length(crv);	# Number of curves
    ncollocpoints=size(collocCoord,1)

    H=complex(zeros(ncollocpoints,ncollocpoints));
    G=complex(zeros(ncollocpoints,ncollocpoints));
    npgauss=12;
    qsi,w=Gauss_Legendre(-1,1,npgauss)
    # qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1 : n #Laço sobre as curvas (elementos)
        uk=unique(crv[i].knots)
        ne=length(uk)-1;
        for j=1:ne
            range=uk[j+1]-uk[j]
            mid=(uk[j+1]+uk[j])/2
            for k=1:ncollocpoints
                xfonte=collocCoord[k,1]
                yfonte=collocCoord[k,2]
                # Integrais de dom�nio
                g,h,id=integra_elem_POT(xfonte,yfonte,crv[i],dcrv[i],range,mid,qsi,w,kmat); # Integra��o sobre o
                #  element (I = int F n.r/r dGama)
                H[k,id+nnos[i]]=H[k,id+nnos[i]]+h;
                G[k,id+nnos[i]]=G[k,id+nnos[i]]+g;
            end
        end
    end
    H=H-E/2;
    return H,G
end

function  integra_elem_POT(xfonte,yfonte,crv,dcrv,range,mid,qsi,w,k)
    # Laplace equation
    # Integra��o sobre os elementos (integral I)
    # N�mero de pontos de Gauss usados na integra��o de I
    npg = length(qsi);
    dudqsi=range/2
    g = complex(zeros(crv.order))
    h = complex(zeros(crv.order))
    id=zeros(Int,crv.order,1)
    # dudqsi=(range[2]-range[1]) / 2;
    # eet = (range[1]+range[2]-2*cpt)/(range[1]-range[2]);
    # eta,Jt=telles(qsi,eet);
    for i = 1 : npg # Percorre os pontos de integra��o
        # qsi_param=convertToParamSpace(eta[i], range);
        p,dp=nrbdeval(crv, dcrv, qsi[i]/2*range+mid)

        B, id = nrbbasisfun(crv,qsi[i]/2*range+mid)
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

        # ZR=real(k*r);
        # Z=complex(0.,ZR);
        # F0C=SpecialFunctions.besselk(0,Z);
        # F1C=SpecialFunctions.besselk(1,Z);
        # qast=-(Z/r*nr*F1C)/(2*pi); #Solução Fundamental da pressão acústica
        # Tast=F0C/(2*pi);    #Solução Fundamental do fluxo de pressão acústica
        #qast = complex(1,0)
        #Tast = complex(1,0)
        h = h + B*qast*dgamadu *dudqsi* w[i]
        g = g + B*Tast*dgamadu *dudqsi* w[i]
    end
    return g,h,id
end

function aplica_CDC(G,H,CDC,E)
    # Aplica as condições de contorno trocando as colunas das matrizes H e G

    ncdc = size(CDC,1); # número de linhas da matriz CDC
    A=copy(H);
    B=copy(G);
    for i=1:ncdc # Laço sobre as condições de contorno
        tipoCDC = CDC[i,2]; # Tipo da condição de contorno
        if tipoCDC == 0 # A temperatura é conhecida
            println(i)
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
    T=complex(zeros(ncdc))
    q=complex(zeros(ncdc))
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

# function calc_phi_pintpot(PONTOS_int,collocCoord,nnos,crv,dcrv,kmat,Tc,qc)
#     n = length(crv);	# Number of curves
#     ncollocpoints = size(collocCoord,1)
#     n_pint=size(PONTOS_int,1)
# 
#     H=complex(zeros(n_pint,ncollocpoints));
#     G=complex(zeros(n_pint,ncollocpoints));
#     npgauss=12;
#     qsi,w=Gauss_Legendre(-1,1,npgauss)
#     # qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
#     for i = 1 : n #Laço sobre as curvas (elementos)
#         uk=unique(crv[i].knots)
#         ne=length(uk)-1;
#         for j=1:ne
#             range=uk[j+1]-uk[j]
#             mid=(uk[j+1]+uk[j])/2
#             for k=1:n_pint
#                 xfonte=PONTOS_int[k,2]
#                 yfonte=PONTOS_int[k,3]
#                 # Integrais de dom�nio
#                 g,h,id=integra_elem_POT(xfonte,yfonte,crv[i],dcrv[i],range,mid,qsi,w,kmat); # Integra��o sobre o
#                 #  element (I = int F n.r/r dGama)
#                 H[k,id+nnos[i]]=H[k,id+nnos[i]]+h;
#                 G[k,id+nnos[i]]=G[k,id+nnos[i]]+g;
#             end
# 
#         end
#     end
# 
#     phi_pint = H*Tc - G*qc;
#     return H,G,phi_pint
# end

function calc_phi_pint(PONTOS_int,collocCoord,nnos,crv,dcrv,kmat,Tc,qc)
    n = length(crv);	# Number of curves
    ncollocpoints = size(collocCoord,1)
    n_pint=size(PONTOS_int,1)

    H=complex(zeros(n_pint,ncollocpoints));
    G=complex(zeros(n_pint,ncollocpoints));
    npgauss=12;
    qsi,w=Gauss_Legendre(-1,1,npgauss)
    # qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1 : n #Laço sobre as curvas (elementos)
        uk=unique(crv[i].knots)
        ne=length(uk)-1;
        for j=1:ne
            range=uk[j+1]-uk[j]
            mid=(uk[j+1]+uk[j])/2
            for k=1:n_pint
                xfonte=PONTOS_int[k,2]
                yfonte=PONTOS_int[k,3]
                # Integrais de dom�nio
                g,h,id=integra_elem(xfonte,yfonte,crv[i],dcrv[i],range,mid,qsi,w,kmat); # Integra��o sobre o
                #  element (I = int F n.r/r dGama)
                H[k,id+nnos[i]]=H[k,id+nnos[i]]+h;
                G[k,id+nnos[i]]=G[k,id+nnos[i]]+g;
            end

        end
    end

    phi_pint = H*Tc - G*qc;
    return phi_pint
end
function calc_phi_pintpot(PONTOS_int,collocCoord,nnos,crv,dcrv,kmat,Tc,qc)
    n = length(crv);	# Number of curves
    ncollocpoints = size(collocCoord,1)
    n_pint=size(PONTOS_int,1)

    H=complex(zeros(n_pint,ncollocpoints));
    G=complex(zeros(n_pint,ncollocpoints));
    npgauss=12;
    qsi,w=Gauss_Legendre(-1,1,npgauss)
    # qsi,w=gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1 : n #Laço sobre as curvas (elementos)
        uk=unique(crv[i].knots)
        ne=length(uk)-1;
        for j=1:ne
            range=uk[j+1]-uk[j]
            mid=(uk[j+1]+uk[j])/2
            for k=1:n_pint
                xfonte=PONTOS_int[k,2]
                yfonte=PONTOS_int[k,3]
                # Integrais de dom�nio
                g,h,id=integra_elem_POT(xfonte,yfonte,crv[i],dcrv[i],range,mid,qsi,w,kmat); # Integra��o sobre o
                #  element (I = int F n.r/r dGama)
                H[k,id+nnos[i]]=H[k,id+nnos[i]]+h;
                G[k,id+nnos[i]]=G[k,id+nnos[i]]+g;
            end

        end
    end

    phi_pint = H*Tc - G*qc;
    return phi_pint
end
function calc_solfundpot(x,y,xd,yd,nx,ny,k)
    # Evaluates the fundamental solutions of the Laplace equation.

    r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
    rx=(x-xd); # x component of the distance
    ry=(y-yd); # y component of the distance
    Tast=-1/(2*pi*k)*log(r); # Fundamental solution for the temperature
    qast=1/(2*pi)*(rx*nx+ry*ny)/r^2; # Fundamental solution for the flux
    return Tast, qast
end



function  calc_solfund(x,y,xd,yd,nx,ny,k)
    # Evaluates the fundamental solutions of the Helmholtz equation.
    r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and field points
    rx=(x-xd)/r; # x component of the distance
    ry=(y-yd)/r; # y component of the distance
    drdn=rx*nx+ry*ny;   # Distance in the normal direction
    ZR=real(k*r);
    Z=complex(0.,ZR);
    F0C=besselk(0,Z);
    F1C=besselk(1,Z);

    qast=-(Z/r*drdn*F1C)/(2*pi); 	# Fundamental solution for the velocity potential
    Tast=F0C/(2*pi);    		# Fundamental solution for the flux
    return Tast,qast
end


