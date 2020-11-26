function Hinterp(indfonte,indbezier,indcoluna,E,Tree1,Tree2,block,crv,kmat,tipoCDC,valorCDC,collocCoord,ninterp=2,compressão=false,ϵ=1e-3)
    # arg = [NOS1,NOS_GEO1,tipoCDC,valorCDC,normal,ELEM1,k]
    #         1      2        3      4        5     6   7
    n = size(block,1)               # Quantidade de Submatrizes
    Aaca = Array{Any}(undef,n,2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    b = complex(zeros(size(valorCDC,1)))       # Cria matriz b, Ax=b, de zeros [Nº de nos x 1]
    for i=1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree1[block[i,1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree2[block[i,2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J
        if block[i,3]==0                # Se esses blocos não são admissiveis
            Aaca[i,1],B = cal_Aeb(indfonte[b1,:],indbezier[b2,:],indcoluna[b2],crv,kmat,tipoCDC,E[b1,:])
            # Salva na linha i da 1º coluna da matriz Aaca a matriz H e salva G na matriz B
            b[b1] = b[b1] + B*valorCDC[vcat(indcoluna[b2]...)]  # Contribuicao para o valor de G*q dos nos que formam b2
        else                              # Caso contrario (Se blocos são admissiveis)
            Aaca[i,1],Aaca[i,2],B=cal_Aeb_interp(indfonte[b1,:],indbezier[b2,:],indcoluna[b2],crv,kmat,tipoCDC,E[b1,:],collocCoord[b1,:],ninterp,compressão,ϵ)
            b[b1] = b[b1] + B*valorCDC[vcat(indcoluna[b2]...)]
        end
        # end

    end
    return Aaca,b
end
function cal_Aeb(indfonte,indbezier,indcoluna, crv, kmat,CDC,E)
    nfonte=size(indfonte,1)
    H = complex(zeros(nfonte,0));
    G = complex(zeros(nfonte,0));
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
        h1 = complex(zeros(nfonte,p+1));
        g1 = complex(zeros(nfonte,p+1));
        for ifonte = 1:nfonte
            k1,k2=indfonte[ifonte,:]
            xfonte = crv[k1].fontes[k2].coords[1]
            yfonte = crv[k1].fontes[k2].coords[2]
            if k1 == i && crv[k1].fontes[k2].pts >= crv[i].range[j,1] && crv[k1].fontes[k2].pts <= crv[i].range[j,2]
                eet = 2 * (crv[k1].fontes[k2].pts - crv[i].range[j,1]) / (crv[i].range[j,2] - crv[i].range[j,1]) - 1
                g, h = integra_sing(xfonte, yfonte, crv[i], qsi, w, crv[i].C[:,:,j], crv[i].conn[j], kmat, eet); # Integra��o sobre o

                if crv[k1].fontes[k2].pts != crv[i].range[j,2]
                    h=h+E[ifonte,indcoluna[ibezier]] /2
                end
            else
                g, h = integra_elem(xfonte, yfonte, crv[i], qsi, w, shapes, derivs, crv[i].C[:,:,j], crv[i].conn[j], kmat); # Integra��o sobre o
            end
            if CDC[i]==1
                h1[ifonte,:] += h;
                g1[ifonte,:] += g;
            else
                h1[ifonte,:] += -g;
                g1[ifonte,:] += -h;
            end

        end
        H= [H h1];
        G= [G g1];
    end
    H , G
end
function cal_Aeb_interp(indfonte,indbezier,indcoluna, crv, kmat,CDC,E,collocCoord,ninterp=2,compressão=true,ϵ=1e-3)
    G = complex(zeros(ninterp*ninterp,0))      # Dimensiona matriz G
    H = complex(zeros(ninterp*ninterp,0))      # Dimensiona matriz H
    xmax=maximum(collocCoord,dims=1)
    xmin=minimum(collocCoord,dims=1)
    xs=criapontosinterp(ninterp)
    n1,n2=calc_fforma(xs)
    xks=n1*xmin+n2*xmax

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
        h1 = complex(zeros(ninterp*ninterp,p+1));
        g1 = complex(zeros(ninterp*ninterp,p+1));
        ci=0
        for i2 =1:ninterp # Laco sobre os pontos fontes
            for i1 =1:ninterp # Laco sobre os pontos fontes
                ci+=1
                xfonte = xks[i1,1]
                yfonte = xks[i2,2]
                g, h = integra_elem(xfonte, yfonte, crv[i], qsi, w, shapes, derivs, crv[i].C[:,:,j], crv[i].conn[j], kmat); # Integra��o sobre o
                if CDC[i]==1
                    h1[ci,:] += h;
                    g1[ci,:] += g;
                else
                    h1[ci,:] += -g;
                    g1[ci,:] += -h;
                end

            end
        end
        H= [H h1];
        G= [G g1];
    end
    fontes=[2*(collocCoord[:,1].-xmin[1])/(xmax[1]-xmin[1]).-1 2*(collocCoord[:,2].-xmin[2])/(xmax[2]-xmin[2]).-1]
    L=lagrange(fontes,xs,ninterp,xs,ninterp)
    if compressão
        H1,H2=recompressão(L,H)
        return H1,H2,L*G
    else
        return L,H,L*G
    end

end




"pg ponto interpolado
x ponto interpolador"
function lagrange(pg,x,n)
    ni = length(pg);
    L = ones(ni,n);
    for j = 1:n
        for i = 1:n
            if (i != j)
                L[:,j] = L[:,j].*(pg .- x[i])/(x[j]-x[i]);
            end
        end
    end
    return L
end

function lagrange(pg,x1,n1,x2,n2)
    l1=lagrange(pg[:,1],x1,n1)
    l2=lagrange(pg[:,2],x2,n2)


    ni=size(pg,1)
    L=zeros(ni,n1*n2)
    for i=1:ni
        L[i,:]=(l1[i,:]*l2[i,:]')[:]
    end
    L
end
function lagrange(pg,x1,n1,x2,n2,x3,n3)
    l1=lagrange(pg[:,1:2],x1,n1,x2,n2)
    l2=lagrange(pg[:,3],x3,n3)


    ni=size(pg,1)
    L=zeros(ni,n1*n2*n3)
    for i=1:ni
        L[i,:]=(l1[i,:]*l2[i,:]')[:]
    end
    L
end
function criapontosinterp(n)
    x = cos.((2*(1:n).-1)*pi/2/n);
end
function calc_fforma(qsi)
    # Calcula as funções de forma lineares contínuas N1 e N2
    N1=1/2*(1.0 .-qsi); # Função de forma N1 => linear contínua
    N2=1/2*(1.0 .+qsi); # Função de forma N2 => linear contínua
    return N1,N2
end



function recompressão(L,H)
    w1,r1=qr(L)
    w1=Array(w1)
    w2,r2=qr(H')
    w2=Array(w2)

    U, S, V = svd(r1*r2')
    ind1=S.>ϵ*S[1]#novo posto
    # ind2=S1.>ϵ*S1[1]#novo posto

    H1=w1[:,ind1]*U[ind1,ind1]*Diagonal(S[ind1])

    H2=V[ind1,ind1]'*w2[:,ind1]'
    H1*H2

    # G1=w1[:,ind2]*U1[ind2,ind2]*diagm(S1[ind2])
    # G2=V1[ind2,ind2]'*w3[:,ind2]'
    # G1*G2
    H1,H2
end
#=
f(x)= 1./(1+25*x.^2);
f(x,y)= 1./(1+25*x.^2+y.^2);
n=30
t =linspace(-1,1,n)
L,x=lagrange(t,11)
xx=t*ones(1,n)
yy=ones(n,1)*t'
ni=10
L2,x2=lagrange2d([xx[:] yy[:]],ni,ni)
A = [ f(x2[i,1],x2[i,2]) for i=1:ni*ni]
scatter(t,t,L2*A)
using Plots
pyplot()
plot(t,f(t))
plot!(t,L*f(x))
surface(xx[:], yy[:],L2*A)
surface(xx[:], yy[:],f)=#

# t1,t2=cal_Aeb(Tree[4],Tree[5],[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
# t3,t4,t5=cal_Aeb_interp(Tree[4],Tree[5],[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
# 100*[(norm(t1)-norm(t3))/norm(t1) (norm(t2)-norm(t4))/norm(t2)]
