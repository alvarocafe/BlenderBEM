# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
#Start----Problem--------------Method-----------------post-processing
#--------------------------------^You are here!----------------------
# Calculate influence matrices with no approximation.
function cal_Aeb(b1,b2,arg)
    # Builds the matrices for the linear system A x = b 
    NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k=arg
    nelem::Int64=size(ELEM)[1]; # Number of elements
    qsitelles,Jtelles = telles(qsi,0); # Evaluating the Telles' points and its Jacobian
    B=complex(zeros(length(b1),length(b2)));	# Allocates matrix B
    A=complex(zeros(length(b1),length(b2)));	# Allocates matrix A
    q=zeros(length(b1),1);	# Allocates array q
    ci=0
    for i in b1 # Loop over the source points
        ci+=1
        xd=NOS[i,2]; # x coordinate of the source point
        yd=NOS[i,3]; # y coordinate of the source point
        cj=0
        for j in b2 # Loop over the elements
            cj+=1
            noi::Int64=ELEM[j,2]; # Initial point of the element
            nof::Int64=ELEM[j,3]; # Final point of the element
            x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
            x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
            y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
            y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
            if i==j # The source point belongs to the element
                g,h = calcula_GeHs(x1,y1,x2,y2,k);	# Singular integration
	        #g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	        #h = 0.5;
                # println("Diferença entre g e gtelles = ", abs(g-gtelles))
                # println("Diferença entre h e htelles = ", abs(h-htelles))
            else # O ponto fonte n�o pertence ao elemento
                g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsi,w,k); # Non singular integration
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
    return A,b
end

function cal_Aebpot(b1,b2,arg)
    # Builds the matrices for the linear system A x = b 
    NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k=arg
    nelem::Int64=size(ELEM)[1]; # Number of elements
    qsitelles,Jtelles = telles(qsi,0); # Evaluating the Telles' points and its Jacobian
    B=complex(zeros(length(b1),length(b2)));	# Allocates matrix B
    A=complex(zeros(length(b1),length(b2)));	# Allocates matrix A
    q=zeros(length(b1),1);	# Allocates array q
    ci=0
    for i in b1 # Loop over the source points
        ci+=1
        xd=NOS[i,2]; # x coordinate of the source point
        yd=NOS[i,3]; # y coordinate of the source point
        cj=0
        for j in b2 # Loop over the elements
            cj+=1
            noi::Int64=ELEM[j,2]; # Initial point of the element
            nof::Int64=ELEM[j,3]; # Final point of the element
            x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
            x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
            y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
            y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
            if i==j # The source point belongs to the element
                #g,h = calcula_GeHs(x1,y1,x2,y2,1,k);	# Singular integration
	        g,h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	        h = -0.5;
                # println("Diferença entre g e gtelles = ", abs(g-gtelles))
                # println("Diferença entre h e htelles = ", abs(h-htelles))
            else
                g,h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsi,w,k); # Non singular integration
            end
	    # Applying the boundary conditions
            if CDC[j,2]==0	# The temperature is known
                B[ci,cj] = -h	# Matrix B receives the value from matrix h
                A[ci,cj] = -g	# Matrix A receives the value from matrix g
            else
                B[ci,cj] = g	# Matrix B receives the value from matrix g
                A[ci,cj] = h	# Matrix A receives the value from matrix h
            end
        end
    end
    b = B*(CDC[b2,3])  # Builds the b array for the linear system
    return A,b
end
function monta_phieq(CDC,x)
    # Apply the boundary conditions on system A x = b, once the vector x has been determined, to separate the values of the velocity potential (phi) and flux (q), described as H phi = G q
    ncdc = length(CDC[:,1]);	# Number of boundary conditions
    nnos = length(x)	# Number of nodes
    phi = complex(zeros(nnos))	# Allocates the vector for the velocity potential
    q = complex(zeros(nnos))	# Allocates the vector for the flux
    for i=1:ncdc # Loop over the boundary conditions
        tipoCDC=CDC[i,2]; # Type of the boundary condition (Neumann or Dirichlet)
        valorCDC=CDC[i,3]; # Value of the boundary condition
        valorcalculado=x[i]; # Value which was previously unknown
        if tipoCDC == 1 # The flux is known
            phi[i] = valorcalculado; # The velocity potential was evaluated
            q[i] = valorCDC; # The flux is given by the boundary condition
        else # The velocity potential is known
            phi[i] = valorCDC; # The velocity potential is given by the boundary coindition
            q[i] = valorcalculado; # The flux was evaluated
        end
    end

    return phi,q
end
function aplica_CDC(G,H,CDC)
    # Applies the boundary conditions by exchanging the columns of matrices H and G and constructs the linear system: Ax = b, where x contains all of the unknowns of the problem.
    ncdc = length(CDC[:,1]); # number of lines in the boundary conditions matrix
    A=H*1.0;	# Right side part of the linear system (unknowns)
    B=G*1.0;	# Left side part of the linear system (knowns)
    for i=1:ncdc # Loop on the number of boundary conditions
        tipoCDC = CDC[i,2]; # Type of the boundary condition
        if tipoCDC == 0 # The velocity potential is known
            colunaA=-A[:,i]; # 
            A[:,i]=-B[:,i]; # Matrix A receives the column from matrix G
            B[:,i]=colunaA; # Matrix b receives the column from matrix H
        end
    end
    valoresconhecidos=CDC[:,3]; # Boundary conditions' values
    b=B*valoresconhecidos; # array b
    return A, b
end


function cal_GeH(NOS,NOS_GEO,ELEM,k,fc,qsi,w)
    # Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.
    
    nelem::Int64=size(ELEM)[1]; # Number of elements
    nnos::Int64=nelem; # Number of nodes
    G=complex(zeros(nnos,nnos)); 	# Allocates matrix G
    H=complex(zeros(nnos,nnos));	# Allocates matrix H
    q=zeros(nnos,1);  # Influence from concentrated sources
    inc = zeros(nnos,1);  # Influence from incident plane waves

    for i=1:nnos # Loop over the source points
        xd=NOS[i,2]; # x coordinate of the source point
        yd=NOS[i,3]; # y coordinate of the source point
        for j=1:nelem # Loop over the elements
            noi::Int64=ELEM[j,2]; # First point of the element
            nof::Int64=ELEM[j,3]; # Second point of the element
            x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
            x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
            y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
            y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
            if i==j # The source point belongs to the element
                g,h = calcula_GeHs(x1,y1,x2,y2,k); 	# Singular integration
	        qsitelles,Jtelles = telles(qsi,-1); # Evaluates the Telles' points and Jacobian
                g1, h1 = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	        qsitelles,Jtelles = telles(qsi,1); # Evaluates the Telles' points and Jacobian
                g2, h1 = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	        qsitelles,Jtelles = telles(qsi,0); # Evaluates the Telles' points and Jacobian
                g4, h4 = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k);	# Singular integration using the Telles transformation
	        h1 = 0.5;
	        g3 = g1 + g2;
                println("Diferença entre g e g3 = ", abs(g-g3))
                println("Diferença entre g e g4 = ", abs(g-g4))
            else # The source point doesnt belong to the element
                g,h = calcula_GeHns(x1,y1,x2,y2,xd,yd,qsi,w,k);	# Non singular integration
            end
            G[i,j] = g
            H[i,j] = h
        end
    end
    # Rigid body consideration for the evaluation of the singular terms of matrix H
    #for m = 1 : nnos
    #    H[m,m] = 0; # let the diagonal be zero
    #    for n = 1 : nnos	# Loop over the source points
    #        if n != m	# If it's not a diagonal term
    #            H[m,m] = H[m,m] - H[m,n]; # Subtracts the term from the diagonal
    #        end;
    #    end;
    #end;

    return G,H
end

function cal_GeHpot(NOS,NOS_GEO,ELEM,k,fc,qsi,w)
    # Evaluates the G and H matrices for the linear system H phi = G q, where phi is a vector containing the values of the velocity potential and q is a vector containing the values of the flux at the boundary of the problem.
    
    nelem::Int64=size(ELEM)[1]; # Number of elements
    nnos::Int64=nelem; # Number of nodes
    G=complex(zeros(nnos,nnos)); 	# Allocates matrix G
    H=complex(zeros(nnos,nnos));	# Allocates matrix H
    q=zeros(nnos,1);  # Influence from concentrated sources
    inc = zeros(nnos,1);  # Influence from incident plane waves
    qsitelles,Jtelles = telles(qsi,0); # Evaluates the Telles' points and Jacobian

    for i=1:nnos # Loop over the source points
        xd=NOS[i,2]; # x coordinate of the source point
        yd=NOS[i,3]; # y coordinate of the source point
        for j=1:nelem # Loop over the elements
            noi::Int64=ELEM[j,2]; # First point of the element
            nof::Int64=ELEM[j,3]; # Second point of the element
            x1=NOS_GEO[noi,2]; # x coordinate of the first point of the element
            x2=NOS_GEO[nof,2]; # x coordinate of the second point of the element
            y1=NOS_GEO[noi,3]; # y coordinate of the first point of the element
            y2=NOS_GEO[nof,3];  # y coordinate of the second point of the element
            if i==j # The source point belongs to the element
                #g,h = calcula_GeHs(x1,y1,x2,y2,k); 	# Singular integration
	        g, h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsitelles,w.*Jtelles,k)
                # Singular integration using the Telles transformation
	        h = -0.5
                #htelles = htelles + 0.5	# Adding the jump term
                #println("Diferença entre g e gtelles = ", abs(g-gtelles))
                #println("Diferença entre h e htelles = ", abs(h-htelles))
            else # The source point doesnt belong to the element
                g,h = calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsi,w,k);	# Non singular integration
            end
            G[i,j] = g
            H[i,j] = h
        end
    end
    # Rigid body consideration for the evaluation of the singular terms of matrix H
    #for m = 1 : nnos
    #    H[m,m] = 0; # let the diagonal be zero
    #    for n = 1 : nnos	# Loop over the source points
    #        if n != m	# If it's not a diagonal term
    #            H[m,m] = H[m,m] - H[m,n]; # Subtracts the term from the diagonal
    #        end;
    #    end;
    #end;

    return G,H
end
function calcula_GeHns(x1,y1,x2,y2,xd,yd,qsi,w,k)
    # Non singular integration
    npg=length(qsi); # Number of integration points
    g=complex(0,0); # Start the sum for g
    h=complex(0,0); # Start the sum for h
    L=sqrt((x2-x1)^2+(y2-y1)^2); # Element length
    dgamadqsi=L/2; # Jacobian
    sx=(x2-x1)/L; # x component of the tangent vector
    sy=(y2-y1)/L; # y component of the tangent vector
    nx=sy; # x component of the normal vector
    ny=-sx; # y component of the normal vector

    for kk=1:npg # Loop over the integration points
        N1,N2 =calc_fforma(qsi[kk]); # Evaluates the shape functions
        x=N1*x1+N2*x2; # Evaluates the x coordinate of the integration point
        y=N1*y1+N2*y2; # Evaluates the y coordinate of the integration point

        Tast,qast =calc_solfund(x,y,xd,yd,nx,ny,k); # Evaluate the fundamental solutions
        #Tast = 1; qast = 1; #To test the integration, let the fundamental solutions be 1. The result of the sum of each line of matrices H and G should be the perimeter of the geometric model.
        h=h+qast*w[kk]*dgamadqsi; # Integration of matrix H
        g=g+Tast*w[kk]*dgamadqsi; # Integration of matrix G
    end
    return g,h
end

function calcula_GeHnspot(x1,y1,x2,y2,xd,yd,qsi,w,k)
    #integração não singular
    n_pint=length(qsi); # Número de pontos de integração (comprimento do
    #    vetor qsi)
    G=0; # Inicializa o somatorio de g
    H=0; # Inicializa o somatorio de h
    for kk=1:n_pint # Laço sobre os pontos de integração
        N1,N2=calc_fforma(qsi[kk]); # Calcula as funções de forma
        dN1dqsi=-0.5
        dN2dqsi=0.5 # Calcula as derivadas das
        #    funções de forma
        x=N1*x1+N2*x2; # Calcula a coordenada x do ponto de integração
        y=N1*y1+N2*y2; # Calcula a coordenada y do ponto de integração

        dxdqsi=dN1dqsi*x1+dN2dqsi*x2;
        dydqsi=dN1dqsi*y1+dN2dqsi*y2;
        dgamadqsi=sqrt(dxdqsi^2+dydqsi^2);

        sx=dxdqsi/dgamadqsi; # Componente x do vetor tangente
        sy=dydqsi/dgamadqsi; # Componente y do vetor tangente
        nx=sy; # Componente x do vetor normal
        ny=-sx; # Componente y do vetor normal

        Tast,qast=calc_solfundpot(x,y,xd,yd,nx,ny,k); # Calcula as soluções fundamentais
        H=H+qast*dgamadqsi*w[kk]; # Integral da matriz H
        G=G+Tast*dgamadqsi*w[kk]; # Integral da matriz G
    end
    return G, H
end

function calcula_GeHs(X1,Y1,X2,Y2,k)
    # integração singular
    #  THIS SUBROUTINE COMPUTES THE VALUES OF THE DIAGONAL
    #  COEFFICIENTS OF THE G MATRIX

    AX=(X2-X1);
    AY=(Y2-Y1);
    SR=sqrt(AX^2+AY^2)/2;
    ZR=real(k*SR);
    Z=complex(0.,ZR);
    S0C=INBESS(Z);
    G=(SR/Z)*S0C/(pi);
    H=0.5;
    return G, H
end

function INBESS(Z)
    #
    #  THIS SUBROUTINE COMPUTES THE INTEGRAL OF THE MODIFIED BESSEL FUNCTION
    #  OF ZERO ORDER BETWEEN 0 AND Z USING A POLINOMIUM.
    #
    #  THE COEFFICIENTS OF THE POLINOMIUM HAVE BEEN OBTAINED FROM THE
    #  TERMS OF THE SERIES EXPANSION

    S0=0.42278434;
    S1=[0.90734120e+01, 0.72756425E+02, 0.25901019E+03,  0.52398212e+03, 0.68598066e+03, 0.62976149e+03,  0.42826262E+03, 0.22451590E+03, 0.93540077E+02,  0.31723213E+02, 0.89292258E+01, 0.21196877E+01,  0.43013488E+00, 0.75474799E-01, 0.11565621E-01,  0.15612002E-02, 0.18705605E-03, 0.20027866E-04,  0.19277885E-05, 0.16772276E-06];

    S2=[0.12000000E+02, 0.64800000E+02, 0.18514286E+03,  0.32400000E+03, 0.38173091E+03, 0.32300308E+03, 0.20566727E+03, 0.10207750E+03, 0.40592223E+02,  0.13221467E+02, 0.35916023E+01, 0.82606852E+00,   0.16293265E+00, 0.27862515E-01, 0.41703893E-02,  0.55091790E-03, 0.64704940E-04, 0.68008195E-05, 0.64341868E-06, 0.55082916E-07];

    R=abs(Z);
    RR=12.0;
    if(R-RR)<0
        N=max(1,min(20.,19.572/log(RR/R)));
    else
        N=20;
        println(" INBESS: RESULTS MAY BE INACCURATE ")
    end
    ZZ=Z/complex(RR,0.);
    ZZ=ZZ*ZZ;
    FINK0=-log(Z/complex(2.,0.))*(complex(1.,0.)+SER(S2,ZZ,1,floor(Int,N)));
    FINK0=Z*(complex(S0,0.)+SER(S1,ZZ,1,floor(Int,N))+FINK0);

    return FINK0
end

function SER(S,ZZ,N1,N2)

    #  THIS FUNCTION COMPUTES THE VALUE OF A SERIES OF N2-N1+1
    #  TERMS CONSISTING OF INCREASING POWERS OF ZZ TIMES THE
    #  COEFFICIENTS IN ARRAY S

    y=complex(0.,0.);
    ZZZ=complex(1.,0.);
    i::Int64 = 1;
    for i=N1:N2
        ZZZ=ZZZ*ZZ;
        y=y+complex(S[i],0.)*ZZZ;
    end

    return y
end
function calc_F(NOS,CW,FR,fc,finc)
    # Evaluates the influence of any concentrated source or incident plane wave of the problem
    nnos=size(NOS,1);	# Number of physical nodes
    q = complex(zeros(nnos));	# Allocates the vector for the concentrated sources
    inc = complex(zeros(nnos));	# Allocates the vector for incident plane waves
    for i = 1:nnos	# Loop over the nodes
        xd=NOS[i,2]; # x coordinate of the source point
        yd=NOS[i,3]; # y coordinate of the source point
        q[i]=calc_q(xd,yd,fc,FR,CW);	# Evaluates the influence of the concentrated sources
        inc[i] = calc_inc(xd,yd,finc,FR,CW);	# Evaluates the influence of the incident plane wave
    end
    return q, inc
end
function calc_inc(xd,yd,finc,k)
    # Evaluates the influence of the incident plane waves described by finc
    # finc = [flag A x y], where flag is either 1 or 0 and indicates the presence of incident plane waves, A is the amplitude of the wave, x and y are the direction of propagation of the wave (sqrt(x^2 + y^2) = 1)
    n_inc = size(finc,1);	# Number of incident waves
    inc = complex(0,0);	# Inlfuence of all incident plane waves
    for i = 1:n_inc	# Loop over the incident waves
        inc += finc[i,4]*exp(complex(0,k*(finc[i,2]*xd + finc[i,3]*yd)));	# Sums the influence of each plane wave
    end
    return inc
end
function calc_q(xd,yd,fc,k)
    # Evaluates the influence from concentrated acoustic sources
    n_inc = size(fc,1);	# Number of acoustic sources
    q = complex(0,0);	# Allocates the influence of the concentrated sources
    for i = 1:n_inc # Loop over the concentrated sources
        x = fc[i,2];  # x coordinate of the source
        y = fc[i,3];  # y coordinate of the source
        # Evaluate the fundamental solution
        r=sqrt((x-xd)^2+(y-yd)^2); # Distance between the source and the field points
        rx=(x-xd)/r; # x component of the distance
        ry=(y-yd)/r; # y component of the distance
        ZR=real(k*r);	# wave number for the distance between the points
        Z=complex(0.,ZR);	# let the wave number be a purely imaginary number
        F0C=SpecialFunctions.besselk(0,Z); 
        Tast=F0C/(2*pi);    # Evaluates the fundamental solution for the flux
        q += fc[i,4]*Tast;	# Evaluates the influence from the concentrated source
    end

    return q
end
function calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,finc,qsi,w,k)
    # Evaluates the velocity potential at internal (or external) points
    n_pint=length(PONTOS_int[:,1]); # Number of internal points
    n_elem=length(phi); # Number of elements
    G_int = complex(zeros(n_pint,n_elem))	# Allocates the matrix for the influence of the flux at the internal points
    H_int = complex(zeros(n_pint,n_elem)) # Allocates the matrix for the influence of the velocity potential at the internal points
    q = complex(zeros(n_pint,1))	# Allocates the array for the influence of concentrated acoustic sources
    inc = complex(zeros(n_pint,1))	# Allocates the array for the influence of incident plane waves

    for i=1:n_pint # Loop over the internal points
        x_fonte=PONTOS_int[i,2]; # x coordinate of the internal point
        y_fonte=PONTOS_int[i,3]; # y coordinate of the internal point
        for j=1:n_elem  # Loop over the elements
            no1::Int64=ELEM[j,2]; # First point of the element
            no2::Int64=ELEM[j,3]; # Second point of the element

            x1=NOS_GEO[no1,2]; # x coordinate of the first point
            y1=NOS_GEO[no1,3]; # y coordinate of the first point

            x2=NOS_GEO[no2,2]; # x coordinate of the second point
            y2=NOS_GEO[no2,3]; # y coordinate of the second point

            G_int[i,j],H_int[i,j] = calcula_GeHns(x1,y1,x2,y2,x_fonte,y_fonte,qsi,w,k) # Non singular integration
        end
        if fc[1,1] > 0 	# If the flag is greater than zero, there are acoustic sources and/or incident plane waves
            q[i,1] = calc_q(x_fonte,y_fonte,fc,k) # Evaluates the influence from acoustic concentrated sources
            inc[i,1] = calc_inc(x_fonte,y_fonte,finc,k) # Evaluates the influence from acoustic incident plane waves
        else
            q[i,1], inc[i,1] =[0 0]	# There is no influence
        end
    end
    #  phi_pint=-(H_int*phi-G_int*qphi-q-inc); # Evaluates the velocity potential for the internal points
    phi_pint=-(H_int*phi-G_int*qphi); # Evaluates the velocity potential for the internal points
    
    return phi_pint
end

function calc_phi_pintpot(PONTOS_int,NOS_GEO,ELEM,phi,qphi,fc,qsi,w,k)
    # Evaluates the velocity potential at internal (or external) points
    n_pint=length(PONTOS_int[:,1]); # Number of internal points
    n_elem=length(phi); # Number of elements
    G_int = complex(zeros(n_pint,n_elem))	# Allocates the matrix for the influence of the flux at the internal points
    H_int = complex(zeros(n_pint,n_elem)) # Allocates the matrix for the influence of the velocity potential at the internal points
    q = complex(zeros(n_pint,1))	# Allocates the array for the influence of concentrated acoustic sources
    inc = complex(zeros(n_pint,1))	# Allocates the array for the influence of incident plane waves

    for i=1:n_pint # Loop over the internal points
        x_fonte=PONTOS_int[i,2]; # x coordinate of the internal point
        y_fonte=PONTOS_int[i,3]; # y coordinate of the internal point
        for j=1:n_elem  # Loop over the elements
            no1::Int64=ELEM[j,2]; # First point of the element
            no2::Int64=ELEM[j,3]; # Second point of the element

            x1=NOS_GEO[no1,2]; # x coordinate of the first point
            y1=NOS_GEO[no1,3]; # y coordinate of the first point

            x2=NOS_GEO[no2,2]; # x coordinate of the second point
            y2=NOS_GEO[no2,3]; # y coordinate of the second point

            G_int[i,j],H_int[i,j] = calcula_GeHnspot(x1,y1,x2,y2,x_fonte,y_fonte,qsi,w,k) # Non singular integration
        end
        if fc[1,1] > 0 	# If the flag is greater than zero, there are acoustic sources and/or incident plane waves
            q[i,1] = calc_q(x_fonte,y_fonte,fc,k) # Evaluates the influence from acoustic concentrated sources
        else
            q[i,1] =0	# There is no influence
        end
    end
    #  phi_pint=-(H_int*phi-G_int*qphi-q-inc); # Evaluates the velocity potential for the internal points
    phi_pint=(H_int*phi-G_int*qphi); # Evaluates the velocity potential for the internal points
    
    return phi_pint
end

function domain_field(xfield,y,T,q,node,dnorm,kmat)
    # -------------------------------------------------------------------------------
    #  Fucntion used to evaluate the temperature/heat flux using the direct
    #     evaluation of the integral representation
    #
    # -------------------------------------------------------------------------------
    if(isempty(xfield))
        f=zeros(1)
        fx=zeros(1)
        fy=zeros(1)
    else
        pi2 = π*2
        nfield=size(xfield,2)
        f=zeros(nfield)
        fx=zeros(nfield)
        fy=zeros(nfield)
        n=size(y,2)
        al= sqrt.((y[1,vec(node[2,:])]-y[1,vec(node[1,:])]).^2+ (y[2,vec(node[2,:])]-y[2,vec(node[1,:])]).^2)

        for j=1:n
            # Comprimento do elemento
            # Dist�ncias entre o ponto interno e os extremos dos elementos
            x11 = y[1,node[1,j]] - xfield[1,:]
            x21 = y[2,node[1,j]] - xfield[2,:]
            x12 = y[1,node[2,j]] - xfield[1,:]
            x22 = y[2,node[2,j]] - xfield[2,:]
            r1 =  sqrt.(x11.^2 + x21.^2); # Dist�ncia para o in�cio do elemento
            r2 =  sqrt.(x12.^2 + x22.^2); # Dist�ncia para o final do elemento

            # Proje��o do vetor dist�ncia no vetor normal ao elemento
            d  =  x11.*dnorm[1,j] + x21.*dnorm[2,j]; # Figura A.1, p�gina 178
            t1 = -x11.*dnorm[2,j] + x21.*dnorm[1,j]; # Dist�ncia T1 da figura A.1
            t2 = -x12.*dnorm[2,j] + x22.*dnorm[1,j]; # Dist�ncia T2 da figura A.1
            ds = abs.(d)
            dtheta = atan2.(ds.*al[j],ds.^2+t1.*t2)

            # Equa��o (A.5) do livro do Liu: elementos da matriz G
            g = -(-dtheta.*ds + al[j] + t1.*real(log.(r1))-t2.*real(log.(r2)))/(pi2*kmat)
            # Equa��o (A.7) com nx=1, ny=0, tx=0, ty=1.
            kkx = (dtheta.*dnorm[1,j] - real(log.(r2 ./r1)).*dnorm[2,j])/(pi2*kmat)
            # Equa��o(A.7) do livro do Liu com nx=0, ny=1, tx=-1, ty=0.
            kky = (dtheta.*dnorm[2,j] + real(log.(r2 ./r1)).*dnorm[1,j])/(pi2*kmat)
            hhx = -(-(t2 ./r2.^2-t1 ./r1.^2).*dnorm[1,j] - d.*(1 ./r2.^2-1 ./r1.^2).*dnorm[2,j])/pi2; # Equa��o (A.8) com nx=1
            # ny=0; tx=0; ty=1.
            hhy = -(-(t2 ./r2.^2-t1 ./r1.^2).*dnorm[2,j] + d.*(1 ./r2.^2-1 ./r1.^2).*dnorm[1,j])/pi2; # Equa��o (A.8) com nx=0
            #     if(d<=0)
            #         dtheta = -dtheta
            #     end
            dtheta = dtheta.*sign.(d)
            h = -dtheta/pi2; # Equa��o (A.6): Elementos da matriz 
            f = f + g*q[j] - h*T[j]; # Integral (2.12) com o termo de
            # dom�nio igual a zero.
            fx = fx + kkx*q[j] - hhx*T[j]; # Integral (2.12) com o termo de
            # dom�nio igual a zero.
            fy = fy + kky*q[j] - hhy*T[j]; # Integral (2.12) com o termo de
            # dom�nio igual a zero.
        end
    end
    return f,fx,fy
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


