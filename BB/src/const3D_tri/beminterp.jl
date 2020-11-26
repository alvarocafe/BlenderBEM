function Hinterp(Tree,block,arg,ninterp=3,compressão=true,ϵ=1e-3)
    # arg = [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k]
    #         1     2     3    4  5  6  7  8
    n = size(block,1)  # Quantidade de Submatrizes
    Aaca = Array{Any}(undef,n,2) # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    b = complex(zeros(size(arg[1],1))) # Cria matriz b, Ax=b, de zeros [Nº de nos x 1]
    for i=1:n # Para cada Submatriz
	# @timeit to "Para cada Submatriz" begin
	b1 = Tree[block[i,1]] # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
	b2 = Tree[block[i,2]] # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
	# Submatriz = Produto cartesiano I x J
	if block[i,3]==0 # Se esses blocos não são admissiveis
	    Aaca[i,1],B = cal_Aeb(b1,b2,arg)
	    # Salva na linha i da 1º coluna da matriz Aaca a matriz H e salva G na matriz B
	    b[b1] = b[b1] + B  # Contribuicao para o valor de G*q dos nos que formam b2
	else  # Caso contrario (Se blocos são admissiveis)
	    Aaca[i,1],Aaca[i,2],L,B=cal_Aeb_interp(b1,b2,arg,ninterp,compressão,ϵ)
	    b[b1] = b[b1] + L*(B*arg[7][b2,3])
	end
	# end
    end
    return Aaca,b
end

function cal_Aeb_interp(b1,b2,arg,ninterp=3,compressão=true,ϵ=1e-3)
    NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k=arg
    nelem::Int64 = size(ELEM)[1]          # Numero de elementos de contorno
    G = complex(zeros(ninterp*ninterp,length(b2)))      # Dimensiona matriz G
    H = complex(zeros(ninterp*ninterp,length(b2)))      # Dimensiona matriz H
    q = complex(zeros(length(b1),1))               # Dimensiona matriz q

    xmax=maximum(NOS[b1,2:4],dims=1)
    xmin=minimum(NOS[b1,2:4],dims=1)
    xs=criapontosinterp(ninterp)
    n1,n2,n3=calc_fforma(xs)
    xks=n1*xmin+n2*xmax+n3*xmin # porque n3*xmin
    ci=0

    for i2 =1:ninterp # Laco sobre os pontos fontes
	for i1 =1:ninterp # Laco sobre os pontos fontes
	    ci+=1
	    xd=xks[i1,1]; # Coordenada x do ponto fonte
	    yd=xks[i2,2]; # Coordenada y do ponto fonte
            zd=xks[i1,3]; # Coordenada z do ponto fonte # por que xks[i1,3]
	    cj=0

	    for j in b2 # Laco sobre os elementos
		cj+=1
		no1::Int64=ELEM[j,2]; # Ponto inicial do elemento
		no2::Int64=ELEM[j,3]; # Ponto final do elemento
		no3::Int64=ELEM[j,4]; # Ponto final do elemento

                x1=NOS_GEO[no1,2]; # Coordenada x do ponto inicial do elemento
		x2=NOS_GEO[no2,2]; # Coordenada x do ponto final do elemento
		x3=NOS_GEO[no3,2]; # Coordenada x do ponto final do elemento
		y1=NOS_GEO[no1,3]; # Coordenada y do ponto inicial do elemento
		y2=NOS_GEO[no2,3];  # Coordenada y do ponto final do elemento
		y3=NOS_GEO[no3,3];  # Coordenada y do ponto final do elemento
		z1=NOS_GEO[no1,4]; # Coordenada z do ponto inicial do elemento
		z2=NOS_GEO[no2,4];  # Coordenada z do ponto final do elemento
		z3=NOS_GEO[no3,4];  # Coordenada z do ponto final do elemento
                g,h = calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k)

		if CDC[j,2]==0
		    G[ci,cj] = -h
		    H[ci,cj] = -g
		else
		    G[ci,cj] = g
		    H[ci,cj] = h
		end
	    end
	end
    end
    if (xmax[1]-xmin[1]) == 0 && (xmax[2]-xmin[2])!=0
	fontes=[0*(NOS[b1,2]-xmin[1])/(xmax[2]-xmin[2])-1 2*(NOS[b1,3]-xmin[2])/(xmax[2]-xmin[2])-1  2*(NOS[b1,4]-xmin[3])/(xmax[3]-xmin[3])-1]
    elseif (xmax[1]-xmin[1]) != 0 && (xmax[2]-xmin[2])!=0
	fontes=[2*(NOS[b1,2]-xmin[1])/(xmax[1]-xmin[1])-1 0*(NOS[b1,3]-xmin[2])/(xmax[1]-xmin[1])-1  2*(NOS[b1,4]-xmin[3])/(xmax[3]-xmin[3])-1]
    elseif (xmax[1]-xmin[1]) != 0 && (xmax[2]-xmin[2])!=0
	fontes=[2*(NOS[b1,2]-xmin[1])/(xmax[1]-xmin[1])-1 2*(NOS[b1,3]-xmin[2])/(xmax[2]-xmin[2])-1  2*(NOS[b1,4]-xmin[3])/(xmax[3]-xmin[3])-1]
    else
	fontes=[0*(NOS[b1,2]-xmin[1])/(1)-1
                0*(NOS[b1,3]-xmin[2])/(1)-1
                2*(NOS[b1,4]-xmin[3])/(xmax[3]-xmin[3])-1]
    end
    L=lagrange(fontes,xs,ninterp,xs,ninterp,xs,ninterp)
    #if b1[1] == 151 && b1[2]==152
    #	println("b1 = $(b1)")
    #	println("fontes = $(fontes)")
    #	println("L = $(L)")
    #	println("NOS[b1,2] = $(NOS[b1,2])")
    #	println("NOS[b1,3] = $(NOS[b1,3])")

    #	println("(xmax[1]-xmin[1]) = $((xmax[1]-xmin[1]))")
    #	println("(xmax[2]-xmin[2]) = $((xmax[2]-xmin[2]))")

    #	println(AAAA)
    #end
    if 1 ==0
	#if compressão
	w1,r1=qr(L)
	w2,r2=qr(H')
	w3,r3=qr(G')

	U, S, V = svd(r1*r2')
	U1, S1, V1 = svd(r1*r3')
	ind1=S.>ϵ*S[1]#novo posto
	# ind2=S1.>ϵ*S1[1]#novo posto
	H1=w1[:,ind1]*U[ind1,ind1]*diagm(S[ind1])
	H2=V[ind1,ind1]'*w2[:,ind1]'
	#H1*H2

	# G1=w1[:,ind2]*U1[ind2,ind2]*diagm(S1[ind2])
	# G2=V1[ind2,ind2]'*w3[:,ind2]'
	# G1*G2
	return H1,H2,L,G
    else
	return L,H,L,G
    end
end


#"pg ponto interpolado
#x ponto interpolador"

function lagrange(pg,x,n)
    ni = length(pg);
    L = ones(ni,n);
    for j = 1:n
	for i = 1:n
	    if (i != j)
		L[:,j] = L[:,j].*(pg - x[i])/(x[j]-x[i]);
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
    x = cos.((2*(1:n) .-1)*pi/2/n);
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
