function cal_Aeb(b1,b2,arg)
  NOS,NOS_GEO,ELEM,xi,w,CDC,k=arg
  nelem::Int64=size(ELEM)[1]; # Número de elementos (n�mero de linhas da matriz ELEM)

  G=zeros(length(b1),length(b2));
  H=zeros(length(b1),length(b2));
  q=zeros(length(b1),1);
  xitri=(xi .+1)./2
  wtri=w./2

  ci=0
  for i in b1 # Laco sobre os pontos fontes
    ci+=1
    xd=NOS[i,2]; 
    yd=NOS[i,3]; 
    zd=NOS[i,4];
    cj=0
    for j in b2
        cj+=1
        nos = ELEM[j,2:4];
        no1=ELEM[j,2]; 
        no2=ELEM[j,3]; 
        no3=ELEM[j,4]; 
   
        x1=NOS_GEO[no1,2]; 
        y1=NOS_GEO[no1,3]; 
        z1=NOS_GEO[no1,4]; 
        
        x2=NOS_GEO[no2,2]; 
        y2=NOS_GEO[no2,3]; 
        z2=NOS_GEO[no2,4]; 
        
        x3=NOS_GEO[no3,2]; 
        y3=NOS_GEO[no3,3]; 
        z3=NOS_GEO[no3,4]; 
        
        n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
      if i==j # O ponto fonte pertence ao elemento
        g,h = calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,xi,w,k)
      else # O ponto fonte n�o pertence ao elemento
        g,h = calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,xitri,wtri,k);
      end
      if CDC[j,2]==0
        G[ci,cj] = -h
        H[ci,cj] = -g
      else
        G[ci,cj] = g
        H[ci,cj] = h
      end
    end

  end
  return H,G
end

function Hinterp(Tree,block,arg,ninterp)
    # arg = [NOS,NOS_GEO,ELEM,qsi,w,CDC,k]
    #         1      2    3   4   5 6   7
    n = size(block,1)               # Quantidade de Submatrizes
    Aaca = Array{Any}(undef,n,2)          # Cria vetor{Any} que armazena submatrizes [Nº de submatrizes x 2]
    b = zeros(size(arg[1],1))       # Cria matriz b, Ax=b, de zeros [Nº de nos x 1]
    t = 0
    t1 = 0
    for i=1:n                       # Para cada Submatriz
        # @timeit to "Para cada Submatriz" begin
        b1 = Tree[block[i,1]]       # Nós I da malha que formam a submatriz (Pontos Fonte) (linhas)
        b2 = Tree[block[i,2]]       # Nós J da malha que formam a submatriz (Pontos Campo) (Colunas)
        # Submatriz = Produto cartesiano I x J
        if block[i,3]==0                # Se esses blocos não são admissiveis
            Aaca[i,1],B = cal_Aeb(b1,b2,arg)
            
            b[b1] = b[b1] + B*arg[6][b2,3] 
        else                              # Caso contrario (Se blocos são admissiveis)
            Aaca[i,1],Aaca[i,2],B=cal_Aeb_interp(b1,b2,arg,ninterp)
 #           A2=Aaca[i,1]*Aaca[i,2]
 #           A,B2 = cal_Aeb(b1,b2,arg)
 #           println("Matriz L")
 #           println(Aaca[i,1])
 #           println("Matriz HH")
 #           println(Aaca[i,2])
 #           println("Matriz exata")
 #           println(A)
 #           println("Matriz interpolada")
 #           println(A2)

            b[b1] = b[b1] + Aaca[i,1]*(B*arg[6][b2,3])

        end
    end
    return Aaca,b
end

function cal_Aeb_interp(b1,b2,arg,ninterp)
    NOS,NOS_GEO,ELEM,xi,w,CDC,k=arg
 #  println("b1=$b1")
 #  println("b2=$b2")
 #  println("xi =$xi")
 #  println("w =$w")
xmax=zeros(1,3)
xmin=zeros(1,3)
xmax[1]=maximum(NOS[b1,2])
xmin[1]=minimum(NOS[b1,2])
xmax[2]=maximum(NOS[b1,3])
xmin[2]=minimum(NOS[b1,3])
xmax[3]=maximum(NOS[b1,4])
xmin[3]=minimum(NOS[b1,4])
npg=length(xi)
ninterp = 5
xs=criapontosinterp(ninterp)
n1,n2=calc_fforma(xs)
xks=n1 .*xmin .+ n2 .*xmax
ϵ=1e-6
xitri=(xi .+1)./2
wtri=w./2

if abs(xmax[1]-xmin[1])<=ϵ && abs(xmax[2]-xmin[2])>ϵ && abs(xmax[3]-xmin[3])>ϵ
    # As coordenadas x são as mesmas, as coordenadas y e z são interpoladas
    fontes=[(2. .*(NOS[b1,3] .-xmin[2])./(xmax[2]-xmin[2]).-1) (2. .*(NOS[b1,4] .-xmin[3])./(xmax[3]-xmin[3]).-1)]
    ninterp3=ninterp # número de pontos usados na interpolação de z
    ninterp2=ninterp # número de pontos usados na interpolação de y
    ninterp1=1 # número de pontos usados na interpolação de x
    L=lagrange(fontes,xs,ninterp,xs,ninterp);
  println("Entrou aqui 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#   println("b1=$b1")
#   println("b2=$b2")
#   println("fontes = $fontes")
#   println("xmin=$xmin")
#   println("xmax=$xmax")
#   println("xs=$xs")
#   println("xks=$xks")
#   println("n1=$n1")
#   println("n2=$n2")
#   println("NOS[b1,1]=$(NOS[b1,1])")
#   println("NOS[b1,2:4]=$(NOS[b1,2:4])")
#   println("NOS[b2,1]=$(NOS[b2,1])")
#   println("NOS[b2,2:4]=$(NOS[b2,2:4])")
#   println("L=$(L[1,:])")
elseif abs(xmax[2]-xmin[2])<=ϵ && abs(xmax[1]-xmin[1])>ϵ && abs(xmax[3]-xmin[3])>ϵ
    # As coordenadas y são as mesmas, as coordenadas x e z são interpoladas
    fontes=[(2. .*(NOS[b1,2] .-xmin[1])./(xmax[1]-xmin[1]).-1) (2. .*(NOS[b1,4] .-xmin[3])./(xmax[3]-xmin[3]).-1)]
    L=lagrange(fontes,xs,ninterp,xs,ninterp);
    ninterp3=ninterp
    ninterp2=1
    ninterp1=ninterp
   println("Entrou aqui 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#   println("b1=$b1")
#   println("b2=$b2")
#   println("fontes = $fontes")
#   println("xmin=$xmin")
#   println("xmax=$xmax")
#   println("xs=$xs")
#   println("xks=$xks")
#   println("n1=$n1")
#   println("n2=$n2")
#   println("NOS[b1,1]=$(NOS[b1,1])")
#   println("NOS[b1,2:4]=$(NOS[b1,2:4])")
#   println("NOS[b2,1]=$(NOS[b2,1])")
#   println("NOS[b2,2:4]=$(NOS[b2,2:4])")
#   println("L=$(L[1,:])")
#   println(" ")

elseif abs(xmax[3]-xmin[3])<=ϵ && abs(xmax[2]-xmin[2])>ϵ && abs(xmax[1]-xmin[1])>ϵ
    # As coordenadas z são as mesmas, as coordenadas x e y são interpoladas
    fontes=[(2. .*(NOS[b1,2] .-xmin[1])./(xmax[1]-xmin[1]).-1) (2. .*(NOS[b1,3] .-xmin[2])./(xmax[2]-xmin[2]).-1)]
    L=lagrange(fontes,xs,ninterp,xs,ninterp);
    ninterp3=1
    ninterp2=ninterp
    ninterp1=ninterp
  println("Entrou aqui 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
elseif abs(xmax[2]-xmin[2])<=ϵ && abs(xmax[3]-xmin[3]) <= ϵ && abs(xmax[1]-xmin[1])>ϵ
    # As coordenadas y e z são as mesmas, as coordenadas x são interpoladas
    fontes=(2. .*(NOS[b1,2] .-xmin[1])./(xmax[1]-xmin[1]).-1);
    ninterp3=1
    ninterp2=1
    ninterp1=ninterp
    L=lagrange(fontes,xs,ninterp);
   println("Entrou aqui 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
elseif abs(xmax[2]-xmin[2])<=ϵ && abs(xmax[1]-xmin[1]) <= ϵ && abs(xmax[3]-xmin[3])>ϵ
    # As coordenadas x e y são as mesmas, as coordenadas z são interpoladas
    fontes=(2. .*(NOS[b1,4] .-xmin[3])./(xmax[3]-xmin[3]).-1);
    L=lagrange(fontes,xs,ninterp);
    ninterp3=ninterp
    ninterp2=1
    ninterp1=1
   println("Entrou aqui 5 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
elseif abs(xmax[1]-xmin[1])<=ϵ && abs(xmax[3]-xmin[3]) <=ϵ && abs(xmax[2]-xmin[2])>ϵ
    # As coordenadas x e z são as mesmas, as coordenadas y são interpoladas
    fontes=(2. .*(NOS[b1,3] .-xmin[2])./(xmax[2]-xmin[2]).-1);
    L=lagrange(fontes,xs,ninterp);
    ninterp3=1
    ninterp2=ninterp
    ninterp1=1
   println("Entrou aqui 6 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
else
    # Todas as coordenadas são interpolaas
    fontes=[(2. .*(NOS[b1,2] .- xmin[1]) ./(xmax[1]-xmin[1]).-1) (2. .*(NOS[b1,3] .-xmin[2])./(xmax[2]-xmin[2]).-1)  (2. .*(NOS[b1,4] .-xmin[3])./(xmax[3]-xmin[3]).-1)];
    ninterp3=ninterp
    ninterp2=ninterp
    ninterp1=ninterp
    L=lagrange(fontes,xs,ninterp,xs,ninterp,xs,ninterp);
#    println(" ")
    println("Entrou aqui 7 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#    println("b1=$b1")
#    println("b2=$b2")
#    println("fontes = $fontes")
#    println("xmin=$xmin")
#    println("xmax=$xmax")
#    println("xs=$xs")
#    println("xks=$xks")
#    println("n1=$n1")
#    println("n2=$n2")
#    println("NOS[b1,1]=$(NOS[b1,1])")
#    println("NOS[b2,1]=$(NOS[b2,1])")
#    println("NOS[b2,2:4]=$(NOS[b2,2:4])")
#    println("L=$(L[1,:])")
end;
    GG = (zeros(ninterp1*ninterp2*ninterp3,length(b2)))      # Dimensiona matriz G
    HH = (zeros(ninterp1*ninterp2*ninterp3,length(b2)))      # Dimensiona matriz H
    G = (zeros(length(b1),length(b2)))      # Dimensiona matriz G
    H = (zeros(length(b1),length(b2)))      # Dimensiona matriz H
    xyz = zeros(ninterp1*ninterp2*ninterp3,3)      # Dimensiona matriz H

ci=0
for i1 =1:ninterp1 # Laco sobre os pontos fontes
    for i2 =1:ninterp2 # Laco sobre os pontos fontes
        for i3 =1:ninterp3 # Laco sobre os pontos fontes
 #       for i3 =1:1 # Laco sobre os pontos fontes
            ci+=1
            xd=xks[i1,1]; # Coordenada x do ponto fonte
            yd=xks[i2,2]; # Coordenada y do ponto fonte
            zd=xks[i3,3]; # Coordenada z do ponto fonte
#	    xd=NOS[b1[ci],2]; 
#	    yd=NOS[b1[ci],3]; 
#	    zd=NOS[b1[ci],4];
            xyz[ci,:]=[xd yd zd]
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
      #          n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3);
      #          g=0
      #          h=0
      #          for l =1:npg
      #              for m =1:npg
      #                  qsi = (1 - xitri[l])*xitri[m]
      #                  N =calc_fformatri(qsi,xitri[l]); #  funções de forma
      #                  x=N[1]*x1+N[2]*x2+N[3]*x3; # coordenada x do ponto de integração
      #                  y=N[1]*y1+N[2]*y2+N[3]*y3 # coordenada y do ponto de integração
      #                  z=N[1]*z1+N[2]*z2+N[3]*z3; # coordenada z do ponto de integração
      #                  J = calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,qsi,xitri[m]);# jacobiano
      #                  Tast,qast=calc_solfund_POT(x,y,z,xd,yd,zd,n,k); # Soluções fundamentais
      #                  g=g+Tast*(1-xitri[l])*w[l]*w[m]*J # Integral da matriz G
      #                  h=h+qast*(1-xitri[l])*w[l]*w[m]*J # Integral da matriz H
      #              end
      #          end
              n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
              g,h = calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,xitri,wtri,k);

                if CDC[j,2]==0
                    GG[ci,cj] = -h
                    HH[ci,cj] = -g
                else
                    GG[ci,cj] = g
                    HH[ci,cj] = h
                end

            end
        end
    end
end

# ci=0
# for i in b1 # Laco sobre os pontos fontes
#   ci+=1
#   xd=NOS[i,2]; 
#   yd=NOS[i,3]; 
#   zd=NOS[i,4];
#   cj=0
#   for j in b2
#       cj+=1
#       nos = ELEM[j,2:4];
#       no1=ELEM[j,2]; 
#       no2=ELEM[j,3]; 
#       no3=ELEM[j,4]; 
#  
#       x1=NOS_GEO[no1,2]; 
#       y1=NOS_GEO[no1,3]; 
#       z1=NOS_GEO[no1,4]; 
#       
#       x2=NOS_GEO[no2,2]; 
#       y2=NOS_GEO[no2,3]; 
#       z2=NOS_GEO[no2,4]; 
#       
#       x3=NOS_GEO[no3,2]; 
#       y3=NOS_GEO[no3,3]; 
#       z3=NOS_GEO[no3,4]; 
#       
#       n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
#       g,h = calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,xitri,wtri,k);
# #     println("CDC[j,2] = $(CDC[j,2]) j = $j")
#      if CDC[j,2]==0
#        G[ci,cj] = -h
#        H[ci,cj] = -g
#      else
#        G[ci,cj] = g
#        H[ci,cj] = h
#      end
#    end

#  end
 
#  println("G")
#  println(G)
 
 
#  println("L*GG")
#  println(L*GG)
 
#  println("H")
#  println(H)

#  println("L*HH")
#  println(L*HH)

#  println("normaH = $(norm(L*HH-H))")
#  println("normaG = $(norm(L*GG-G))")
#  println(" ")
#  println(" ")
#  println(" ")


   
 #  if (norm(L*xyz - NOS[b1,2:4 ]) > 1e-8)
 #     println(" ")
 #     println("Erro de interpolação !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
 #     println("L*xyz=$(L*xyz)")
 #     println("NOS[b1,1]=$(NOS[b1,1])")
 #     println("NOS[b1,2:4]=$(NOS[b1,2:4])")
 #     println("xyz=$xyz")
 #     println("norm = $(norm(L*xyz-NOS[b1,2:4]))")
 #     println(" ")
  # end
    return L,HH,GG
end

function criapontosinterp(n)
    x= cos.((2. .*(n:-1:1) .-1).*pi./2. ./n)
end

function lagrange(pg,x,n)
    ni = length(pg);
    L = ones(ni,n);
    for k=1:ni
    for j = 1:n
        for i = 1:n
            if (i != j)
                L[k,j] = L[k,j]*(pg[k] - x[i])/(x[j]-x[i]);
            end
        end
    end
    end
    return L
end


function lagrange(pg,x1,n1,x2,n2)
    l2=lagrange(pg[:,1],x1,n1)
    l3=lagrange(pg[:,2],x2,n2)
    ni=size(pg,1)
    L=zeros(ni,n1*n2)
    for i=1:ni
        ii=0
 #       for j=1:ninterp
           for k=1:n1
            for l=1:n2
                ii=ii+1
#                L[i,ii]=l1[i,j]*l2[i,k]*l3[i,l]
                L[i,ii]=l2[i,k]*l3[i,l]
            end
           end
#         end
    end
    return L
end


#function lagrange(pg,x1,n1,x2,n2,x3,n3)                                   
#    l1=lagrange(pg[:,1:2],x1,n1,x2,n2)                                    
#    l2=lagrange(pg[:,3],x3,n3)                                                                          
#    ni=size(pg,1)                                                         
#    L=zeros(ni,n1*n2*n3)                                                  
#    for i=1:ni                                                            
#        L[i,:]=(l1[i,:]*l2[i,:]')[:]                                      
#    end                                                                   
#    L                                                                     
#end 


function lagrange(pg,x1,n1,x2,n2,x3,n3)                                   
    l1=lagrange(pg[:,1],x1,n1)
    l2=lagrange(pg[:,2],x2,n2)
    l3=lagrange(pg[:,3],x3,n3)
    ni=size(pg,1)
    L=zeros(ni,n1*n2*n3)
    for i=1:ni
        ii=0
        for j=1:n1
           for k=1:n2
            for l=1:n3
                ii=ii+1
                L[i,ii]=l1[i,j]*l2[i,k]*l3[i,l]
            end
           end
         end
    end
    return L
end





function matvec(hmat,b,block,Tree)
  v=b*0
  for i =1:length(block[:,3])
    if block[i,3]==1
      v[Tree[block[i,1]]]+=hmat[i,1]*(hmat[i,2]*b[Tree[block[i,2]]])
    else
      v[Tree[block[i,1]]]+=hmat[i,1]*b[Tree[block[i,2]]]
    end
  end
  v
end

function montacheia(hmat,block,Tree,n)
A=zeros(n,n)
  for i =1:length(block[:,3])
    if block[i,3]==1
      A[Tree[block[i,1]],Tree[block[i,2]]]=hmat[i,1]*hmat[i,2]
    else
      A[Tree[block[i,1]],Tree[block[i,2]]]=hmat[i,1]
    end
  end
  A
end
