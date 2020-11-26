function cal_Aeb(b1,b2,arg)
# Builds the matrices for the linear system A x = b
  NOS,NOS_GEO,ELEM,k,qsi,w,inc,CDC = arg
  #NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k=arg
  nelem::Int64=size(ELEM)[1]; # Number of elements
  qsitelles,Jtelles = telles(qsi,0); # Evaluating the Telles' points and its Jacobian
  npg=4; # Number of integration points
  qsiquad,wquad = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
  B=complex(zeros(length(b1),length(b2)));	# Allocates matrix B
  A=complex(zeros(length(b1),length(b2)));	# Allocates matrix A
  q=zeros(length(b1),1);	# Allocates array q
  ci=0
  for i in b1 # Loop over the source points
    ci+=1
    xd=NOS[i,2]; # x coordinate of the source point
    yd=NOS[i,3]; # y coordinate of the source point
    zd=NOS[i,3]; # y coordinate of the source point

    cj=0
    for j in b2 # Loop over the elements
      cj+=1
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

      n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); # vetor unit�rio normal ao elemento
      if i==j # The source point belongs to the element
        #g,h = calcula_GeHs(x1,y1,x2,y2,1,k);	# Singular integration
	#g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	#h = 0.5;
    G[i,j],H[i,j]= calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,k)#
    g,h= calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsiquad,wquad,k)# Integra��o singular

        # println("Diferença entre g e gtelles = ", abs(g-gtelles))
        # println("Diferença entre h e htelles = ", abs(h-htelles))
      else # O ponto fonte n�o pertence ao elemento
	    G[i,j],H[i,j]=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o singular
        g,h=calcula_GeHns(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi,w,k); # Integra��o singular
        #g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsi,w,k); # Non singular integration
      end
	# Applying the boundary conditions
      if CDC[j,2]==0	# The velocity potential is known
        B[ci,cj] = -h	# Matrix B receives the value from matrix h
        A[ci,cj] = -g	# Matrix A receives the value from matrix g
      else
        B[ci,cj] = g	# Matrix B receives the value from matrix g
        A[ci,cj] = h	# Matrix A receives the value from matrix h
      end
    end
  end
	b = B*(CDC[b2,3])  # Builds the b array for the linear system
return A,b,G,H
end
