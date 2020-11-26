function CalcHeG(nnos, crv, kmat)
    dcrvs = map(x -> nrbderiv(x), crv)
    n = length(crv);	# Number of curves
    ncollocpoints = size(collocCoord, 1)

    E = zeros(ncollocpoints, ncollocpoints);
    H = zeros(ncollocpoints, ncollocpoints);
    G = zeros(ncollocpoints, ncollocpoints);
    npgauss = 12;
    qsi, w = gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for i = 1:n
        p = crv[i].order - 1
        shapes  = zeros(npgauss, (p + 1));
        derivs  = zeros(npgauss, (p + 1), 2);
        for gp = 1:size(w, 1)
            shapes[gp,:], derivs[gp,:,:] = bernsteinbasis(p, 0, qsi[gp], 0);
        end
        for j = 1:size(crv[i].conn, 1)
            k = 1
            for k1 = 1:n
                for k2 = 1:crv[k1].number
                    xfonte = crv[k1].fontes[k2].coords[1]
                    yfonte = crv[k1].fontes[k2].coords[2]
                    if k1 == i && crv[k1].fontes[k2].pts >= crv[i].range[j,1] && crv[k1].fontes[k2].pts <= crv[i].range[j,2]
                        eet = 2 * (crv[k1].fontes[k2].pts - crv[i].range[j,1]) / (crv[i].range[j,2] - crv[i].range[j,1]) - 1
                        g, h = integra_sing(xfonte, yfonte, crv[i], qsi, w, crv[i].C[:,:,j], crv[i].conn[j], kmat, eet); # Integra��o sobre o
                        h-=crv[k1].fontes[k2].basis/2
                    else
                        g, h = integra_elem(xfonte, yfonte, crv[i], qsi, w, shapes, derivs, crv[i].C[:,:,j], crv[i].conn[j], kmat); # Integra��o sobre o
                    end
                    H[k,crv[i].conn[j] .+ nnos[i]] += h;
                    G[k,crv[i].conn[j] .+ nnos[i]] += g;
                    k += 1
                end
            end
        end
    end
    k = 1
    for k1 = 1:n
        uu = unique(crv[k1].knots)
        for k2 = 1:crv[k1].number
            ind = sum(crv[k1].fontes[k2].pts .> uu)
            E[k,crv[k1].conn[ind] .+ nnos[k1]] += crv[k1].fontes[k2].basis
            k += 1
        end
    end
    H, G, E
end
function  integra_elem(xfonte, yfonte, crv, qsi, w, shapes, derivs, C, conn, k)
    # Integra��o sobre os elementos (integral I)
    g = zeros(crv.order)
    h = zeros(crv.order)
    for i = 1:size(w, 1) # Percorre os pontos de integra��o
        R, dRdxi = basisfundecomp(shapes[i,:], derivs[i,:,:], C, crv.coefs[4,conn])
        # derivadas das fun��es de forma
        cf = crv.coefs[1:2,conn] ./ [crv.coefs[4,conn]';crv.coefs[4,conn]']
        p = cf * R
        dp = cf * dRdxi;
        dgamadqsi = norm(dp)
        nx = dp[2] / dgamadqsi; # Componente x do vetor normal unit�rio
        ny = -dp[1] / dgamadqsi; # Componente y do vetor normal unit�rio
        x = p[1] - xfonte # Calcula a coordenada x do ponto de integra��o
        y = p[2] - yfonte # Calcula a coordenada y do ponto de integra��o
        r = √(x^2 + y^2); # raio
        rx = x / r; # Componente x do vetor unit�rio r
        ry = y / r; # Componente y do vetor unit�rio r
        nr = nx * rx + ny * ry; # Produto escalar dos vetores unit�rios r e n
        Tast = -1 / (2 * pi * k) * log(r)
        qast = 1 / (2 * pi) * (rx * nx + ry * ny) / r

        h = h + R * qast * dgamadqsi * w[i]
        g = g + R * Tast * dgamadqsi * w[i]
    end
    return g, h
end

function  integra_sing(xfonte, yfonte, crv, qsi, w, C, conn, k, eet)
    # Integra��o sobre os elementos (integral I)
    eta, Jt = telles(qsi, eet)
    g = zeros(crv.order)
    h = zeros(crv.order)
    for i = 1:size(w, 1) # Percorre os pontos de integra��o
        shapes, derivs = bernsteinbasis(crv.order - 1, 0, eta[i], 0);
        R, dRdxi = basisfundecomp(shapes, derivs, C, crv.coefs[4,conn])
        # derivadas das fun��es de forma
        cf = crv.coefs[1:2,conn] ./ [crv.coefs[4,conn]';crv.coefs[4,conn]']
        p = cf * R
        dp = cf * dRdxi;
        dgamadqsi = norm(dp)
        nx = dp[2] / dgamadqsi; # Componente x do vetor normal unit�rio
        ny = -dp[1] / dgamadqsi; # Componente y do vetor normal unit�rio
        x = p[1] - xfonte # Calcula a coordenada x do ponto de integra��o
        y = p[2] - yfonte # Calcula a coordenada y do ponto de integra��o
        r = √(x^2 + y^2); # raio
        rx = x / r; # Componente x do vetor unit�rio r
        ry = y / r; # Componente y do vetor unit�rio r
        nr = nx * rx + ny * ry; # Produto escalar dos vetores unit�rios r e n
        Tast = -1 / (2 * pi * k) * log(r)
        qast = 1 / (2 * pi) * (rx * nx + ry * ny) / r
        h = h + R * qast * dgamadqsi * w[i] * Jt[i]
        g = g + R * Tast * dgamadqsi * w[i] * Jt[i]
    end
    return g, h
end
function aplica_CDC(G, H, CDC, E)
    # Aplica as condições de contorno trocando as colunas das matrizes H e G
    ncdc = size(CDC, 1); # número de linhas da matriz CDC
    A = H;
    B = G;
    for i = 1:ncdc # Laço sobre as condições de contorno
        tipoCDC = CDC[i,2]; # Tipo da condição de contorno
        if tipoCDC == 0 # A temperatura é conhecida
            colunaA = -A[:,i]; # Coluna da matriz H que será trocada
            A[:,i] = -B[:,i]; # A matriz H recebe a coluna da matriz G
            B[:,i] = colunaA; # A mstriz G recebe a coluna da matriz H
        end
    end
    valoresconhecidos = E \ CDC[:,3] # Valores das condições de contorno
    b = B * valoresconhecidos; # vetor b
    return A, b
end

function monta_Teq(tipoCDC,valorCDC,crv, x)
    # Separa fluxo e temperatura
    # ncdc = número de linhas da matriz CDC
    # T = vetor que contém as temperaturas nos nós
    # q = vetor que contém o fluxo nos nós
    n = length(crv);    # Number of curves
    T=0*x
    q=0*x
    z=1
    for k=1:n
        for i=1:crv[k].number
            if tipoCDC[k] == 0
                T[z]=valorCDC[z]
                q[z]=x[z]
            else
                T[z]=x[z]
                q[z]=valorCDC[z]
            end
            z+=1
        end
    end
    return T, q
end

function monta_Teq(tipoCDC,valorCDC, x)
    # Separa fluxo e temperatura
    # ncdc = número de linhas da matriz CDC
    # T = vetor que contém as temperaturas nos nós
    # q = vetor que contém o fluxo nos nós
    T = 0*x
    q = 0*x
    for ibezier = 1:size(indbezier,1) # Laço sobre as condições de contorno
        valorcalculado = x[i] # Valor que antes era desconhecido
        if tipoCDC == 1 # Fluxo é conhecido
            T[i] = valorcalculado; # A temperatura é o valor calculado
            q[i] = valorCDC; # O fluxo é a condiçao de contorno
        else # A temperatura é conhecida
            T[i] = valorCDC; # A temperatura é a condiçao de contorno
            q[i] = valorcalculado; # O fluxo é o valor calculado
        end
    end
    return T, q
end


function indices(crv)
    n = length(crv);    # Number of curves

    z=0;#ncollocpoints
    for k=1:n
        for i=1:crv[k].number
            z=z+1
        end
    end

    numcurva=zeros(Integer,z)
    collocPts=zeros(z)
    tCDC=zeros(0,1)
    collocCoord=zeros(z,2)
    nnos=zeros(Integer,n)
    for k=1:n
        p=crv[k].order-1;
        nnos[k]=crv[k].number;
        valorCDC=CCSeg[k,3];
        tipoCDC=CCSeg[k,2];
        tCDC = [tCDC;tipoCDC];
    end
    vCDC=zeros(0,1)

    for k=1:n
        p=crv[k].order-1;
        nnos[k]=crv[k].number;
        valorCDC=CCSeg[k,3];
        tipoCDC=CCSeg[k,2];
        for i=1:crv[k].number
            vCDC= [vCDC;valorCDC];
        end
    end

    nnos2=cumsum([0 nnos'],dims=2);


    indfonte = Array{Int}(undef,z, 2);
    indbezier = Array{Int}(undef,0, 2);
    indcoluna = Array{Any}(undef,0,1)


    for i = 1:n
        for j = 1:size(crv[i].conn, 1)
            indbezier=[indbezier
                       i j]
            indcoluna=[indcoluna
                       [(crv[i].conn[j] .+ nnos2[i])]]
        end
    end
    k = 1
    for k1 = 1:n
        for k2 = 1:crv[k1].number
            indfonte[k,:]=[k1 k2]
            k += 1
        end
    end




    E = zeros(z,z);

    k = 1
    for k1 = 1:n
        uu = unique(crv[k1].knots)
        for k2 = 1:crv[k1].number
            ind = sum(crv[k1].fontes[k2].pts .> uu)
            E[k,crv[k1].conn[ind] .+ nnos2[k1]] += crv[k1].fontes[k2].basis
            k += 1
        end
    end
    cont=1
    for c in crv
        for f in c.fontes
            collocCoord[cont,:]=f.coords[1:2]
            collocPts[cont]=f.pts
            cont+=1
        end
    end


    indfonte,indcoluna,indbezier,tCDC,vCDC,E,collocCoord,collocPts
end



function CalcAeb(indfonte,indcoluna,indbezier, crv, kmat,E,CDC)
    n = length(crv);    # Number of curves
    ncollocpoints = size(collocCoord, 1)

    H = zeros(ncollocpoints, ncollocpoints);
    G = zeros(ncollocpoints, ncollocpoints);
    npgauss = 36;
    qsi, w = gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for ibezier = 1:size(indbezier,1)
        i,j=indbezier[ibezier,:]
        p = crv[i].order - 1
        shapes  = zeros(npgauss, (p + 1));
        derivs  = zeros(npgauss, (p + 1), 2);
        for gp = 1:size(w, 1)
            shapes[gp,:], derivs[gp,:,:] = bernsteinbasis(p, 0, qsi[gp], 0);

        end
        for ifonte = 1:size(indfonte,1)
            k1,k2=indfonte[ifonte,:]
            xfonte = crv[k1].fontes[k2].coords[1]
            yfonte = crv[k1].fontes[k2].coords[2]
            if k1 == i && crv[k1].fontes[k2].pts >= crv[i].range[j,1] && crv[k1].fontes[k2].pts <= crv[i].range[j,2]
                eet = 2 * (crv[k1].fontes[k2].pts - crv[i].range[j,1]) / (crv[i].range[j,2] - crv[i].range[j,1]) - 1
                g, h = integra_sing(xfonte, yfonte, crv[i], qsi, w, crv[i].C[:,:,j], crv[i].conn[j], kmat, eet); # Integra��o sobre o
                # h=h-crv[k1].fontes[k2].basis/2
                if crv[k1].fontes[k2].pts != crv[i].range[j,2]
                    h=h-E[ifonte,indcoluna[ibezier]] /2
                end
                # @show h
                # @show crv[k1].fontes[k2].basis
            else
                g, h = integra_elem(xfonte, yfonte, crv[i], qsi, w, shapes, derivs, crv[i].C[:,:,j], crv[i].conn[j], kmat); # Integra��o sobre o
                # h=0*h#-crv[k1].fontes[k2].basis/2
            end
            if CDC[i]==1
                H[ifonte,indcoluna[ibezier]] += h;
                G[ifonte,indcoluna[ibezier]] += g;
            else
                H[ifonte,indcoluna[ibezier]] += -g;
                G[ifonte,indcoluna[ibezier]] += -h;
            end
        end
   end
    return H , G
end

function calc_pintpot(PONTOS_int,indcoluna,indbezier, crv, kmat,desl,tra)
    n = length(crv);    # Number of curves
    ncollocpoints = size(collocCoord, 1)
    n_p_int=size(PONTOS_int,1)
    H = complex(zeros(n_p_int, ncollocpoints));
    G = complex(zeros(n_p_int, ncollocpoints));
    npgauss = 12;
    qsi, w = gausslegendre(npgauss) # Calcula pesos e pontos de Gauss
    for ibezier = 1:size(indbezier,1)
        i,j=indbezier[ibezier,:]
        p = crv[i].order - 1
        shapes  = zeros(npgauss, (p + 1));
        derivs  = zeros(npgauss, (p + 1), 2);
        for gp = 1:size(w, 1)
            shapes[gp,:], derivs[gp,:,:] = bernsteinbasis(p, 0, qsi[gp], 0);

        end
        for ifonte = 1:size(PONTOS_int,1)
            xfonte = PONTOS_int[ifonte,2]
            yfonte = PONTOS_int[ifonte,3]
            g, h = integra_elem(xfonte, yfonte, crv[i], qsi, w, shapes, derivs, crv[i].C[:,:,j], crv[i].conn[j], kmat); # Integra��o sobre o
            H[ifonte,indcoluna[ibezier]] += h;
            G[ifonte,indcoluna[ibezier]] += g;   
        end
     
    end
    desl_pint=H*desl-G*tra
    return desl_pint
end
