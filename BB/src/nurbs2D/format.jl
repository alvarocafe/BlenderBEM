function format_dad(PONTOS,SEGMENTOS,MALHA,BCSeg)
    # Programa para formata��o dos dados de entrada
    # [NOS_GEO,NOS,ELEM,tipoCDC,valorCDC,normal]
    #   Autor: Lucas Campos
    #include("dad_0.jl")
    num_elementos=sum(MALHA[:,2])
    NOS_GEO=zeros(num_elementos,2)
    ELEM=zeros(UInt32,num_elementos,2)
    NOS=zeros(num_elementos,2)
    valorCDC=zeros(num_elementos)
    tipoCDC=falses(num_elementos)
    normal=zeros(2,num_elementos)
    cont_nos = 0;  # Counter to the physical nodes
    cont_el = 0;	# Counter to the elements (the number of physical and geometric elements is the same).
    num_lin = length(SEGMENTOS[:,1]);	# N�mero de linhas no contorno
    p_ini = SEGMENTOS[1,2];
    crv=Array{Curve}(num_lin)

    #______________________________________________________________________
    # Definition of the biggest dimension of the problem
    max_dl = 0
    for lin = 1 : num_lin
        p1 = convert(Int32,SEGMENTOS[lin,2]);
        p2 = convert(Int32,SEGMENTOS[lin,3]);
        xp1 = PONTOS[p1,2]
        yp1 = PONTOS[p1,3]
        xp2 = PONTOS[p2,2]
        yp2 = PONTOS[p2,3]
        dl = √((xp1-xp2)^2+(yp1-yp2)^2)
        if dl > max_dl
            max_dl = dl
        end
    end
    #_____________________________________________________________________

    no_ini=1
    t=1
    p2=0
    no1_prox=0
    while(t<=num_lin)  	# While over all lines
        # Coordinates of the initial and final PONTOS of each line
        # [x1l,y1l,x2l,y2l)
        p1  = convert(Int32,SEGMENTOS[t,2])
        p2  = convert(Int32,SEGMENTOS[t,3])
        x1l = PONTOS[p1,2]
        y1l = PONTOS[p1,3]
        x2l = PONTOS[p2,2]
        y2l = PONTOS[p2,3]
        if(SEGMENTOS[t,4]==0) # The segment is a straight line
            # Increment in x and y direction
            delta_x = x2l - x1l
            delta_y = y2l - y1l
            crv[t]= nrbmak([x1l x2l;y1l y2l; 0 0], [0, 0, 1,1]);
        else #The segment is an arc
            # Compute the center of the arc and its coordinates
            r = SEGMENTOS[t,4]
            xc,yc=calcula_centroiso(x1l,y1l,x2l,y2l,r)
            # Distance between p1 and c (r1) and between p2 and c (r2)
            r1 = √((x1l-xc)^2+(y1l-yc)^2)
            r2 = √((x2l-xc)^2+(y2l-yc)^2)
            if abs(r1-r2)<.00001*max_dl
                # Compute the angle between the lines from point c to p1 [tet1) and c to p2 (tet2]
                tet1,tet2 = calcula_arcoiso(x1l,y1l,x2l,y2l,xc,yc)
                if tet2 < tet1
                    tet2 = tet2 + 2*π
                end

                # Angle of the sector defined by the arc
                if SEGMENTOS[t,4] > 0
                    tet = abs(tet2-tet1)
                    sig = 1.
                else
                    tet = 2*π-abs(tet2-tet1)
                    sig = -1.
                end
                delta=1e-12
                if abs(tet)-delta <= pi/2
                    narcs = 1;                # number of arc segments
                    knots = [0, 0, 0, 1, 1, 1];
                elseif abs(tet)-delta <= pi
                    narcs = 2;
                    knots = [0, 0, 0, 0.5, 0.5, 1, 1, 1];
                elseif abs(tet)-delta <= 3*pi/2
                    narcs = 3;
                    knots = [0, 0, 0, 1/3, 1/3, 2/3, 2/3, 1, 1, 1];
                else
                    narcs = 4;
                    knots = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1];
                end

                dtet = tet/(2*narcs);     #arc segment tet angle/2

                # determine middle control point and weight
                wm = cos(dtet);
                x  = r*wm;
                y  = r*sin(dtet);
                xm = x+y*tan(dtet);

                # arc segment control points
                ctrlpt = [ x wm*xm x;    # w*x - coordinate
                           -y 0     y;    # w*y - coordinate
                           1 wm    1];   # w   - coordinate

                # build up complete arc from rotated segments
                coefs = zeros(3,2*narcs+1);   # nurbs control points of arc
                angle=tet1 + dtet
                sn = sin(angle);
                cn = cos(angle);
                coefs[:,1:3] = [cn -sn 0 ; sn cn 0 ; 0 0 1]*ctrlpt;     # rotate to start angle
                angle=2*dtet
                sn = sin(angle);
                cn = cos(angle);
                for n = 2:narcs
                    m = 2*n+[0, 1];
                    coefs[:,m] = [cn -sn 0 ; sn cn 0 ;0 0 1]*coefs[:,m-2];
                end

                # vectrans arc if necessary
                coefs = [1 0 xc; 0 1 yc; 0 0 1]*coefs;
		linha,coluna=size(coefs)
		coefsz=[coefs[1:2,:];zeros(coluna)';coefs[3,:]']
                crv[t] = nrbmak(coefsz,knots)
            else
                error("Error in the data input file: Wrong central point")
            end
        end
        t=t+1
    end                                  # end of while(t<num_lin)
        dcrv=map(x->nurbs2D.nrbderiv(x),crv)
    n = length(crv)
    z=0;
    for k=1:n
        for i=1:crv[k].number
            z=z+1
        end
    end
    numcurva=zeros(Integer,z)
    collocPts=zeros(z)
    CDC=zeros(z,3)
    collocCoord=zeros(z,3)
    z=0;
    nnos=zeros(Integer,n)
    for k=1:n
        p=crv[k].order-1;
        nnos[k]=crv[k].number;
        valorCDC=BCSeg[k,3];
        tipoCDC=BCSeg[k,2];
        for i=1:crv[k].number
            z=z+1;
            numcurva[z]=k;
            collocPts[z]=sum(crv[k].knots[(i+1):(i+p)])/p;
            if(i==2)
                collocPts[z-1]=(collocPts[z]+collocPts[z-1])/2;
            end
            if(i==nnos[k])
                collocPts[z]=(collocPts[z]+collocPts[z-1])/2;
            end

            CDC[z,:] = [z,tipoCDC,valorCDC];
        end
    end
    nnos2=cumsum([0 nnos'],2)
    E=zeros(length(collocPts),length(collocPts));
    for i=1:length(collocPts)
        collocCoord[i,:]=nurbs2D.nrbeval(crv[numcurva[i]], collocPts[i]);
        B, id = nurbs2D.nrbbasisfun(crv[numcurva[i]],collocPts[i])
        E[i,id+nnos2[numcurva[i]]]=B
    end
    return collocCoord,nnos2,crv,dcrv,E,CDC
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

function calcula_arcoiso(x1,y1,x2,y2,xc,yc)
    # Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
    # horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
    # and the horizontal direction


    dx1 = x1 - xc; dy1 = y1 - yc
    dx2 = x2 - xc; dy2 = y2 - yc

    # Computation of tet1
    if dy1 == 0				# The point 1 and the center have the same y coordinate
        if x1 > xc
            tet1 = 0
        else  # (x1 < xc)
            tet1 = π
        end
    elseif dx1 == 0				# The point 1 and the center have the same x coordinate
        if y1 > yc
            tet1 = π/2
        else  # (y1 < yc)
            tet1 = -π/2
        end
    else  # (dx1~=0 e dy1~=0)
        tet1 = atan(dy1/dx1);
        if dx1<0 && tet1<0
            tet1 = π + tet1
        elseif dx1 < 0 && tet1>0
            tet1 = -π + tet1
        end
    end

    # Computation of tet2
    if dy2 == 0				# The point 2 and the center have the same y coordinate
        if x2 > xc
            tet2 = 0
        else  # (x2 < xc)
            tet2 = π
        end
    elseif dx2 == 0				# The point 2 and the center have the same x coordinate
        if y2 > yc
            tet2 = π/2
        else  # (y2 < yc)
            tet2 = -π/2
        end
    else  # (dx2~=0 e dy2~=0)
        tet2 = atan(dy2/dx2);
        if dx2<0 && tet2<0
            tet2 = π + tet2
        elseif dx2 < 0 && tet2>0
            tet2 = -π + tet2
        end
    end
    [tet1,tet2]
end


function calcula_centroiso(x1,y1,x2,y2,raio)
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

function monta_Teqiso(tipoCDC,valorCDC,x)
    # Separa fluxo e temperatura
    
    # ncdc = n�mero de linhas da matriz CDC
    # T = vetor que cont�m as temperaturas nos n�s
    # q = vetor que cont�m o fluxo nos n�s
    
    ncdc = length(tipoCDC)
    T=zeros(ncdc)
    q=zeros(ncdc)
    for i=1:ncdc # La�o sobre as condi��es de contorno
        if tipoCDC[i] == 1 # Fluxo � conhecido
            T[i] = x[i]; # A temperatura � o valor calculado
            q[i] = valorCDC[i]; # O fluxo � a condi�ao de contorno
        else # A temperatura � conhecida
            T[i] = valorCDC[i]; # A temperatura � a condi�ao de contorno
            q[i] = x[i]; # O fluxo � o valor calculado
        end
    end
    return T,q
end
function mostra_geo(crvs)
    p=plot(legend=:none,aspect_ratio=:equal)
    for i in crvs
        p=nrbplot(i)
    end
    p
end
#NURBS toolbox

type Curve
    knots
    coefs
    number
    order
end

function convertToParamSpace( xi, range)
    # Take in a coordinate in the parent coordinate system xi in [-1,1]
    # and convert it to the parameter space
    
    xi_param=( (range[2] - range[1]) * xi  + (range[2] + range[1] ) ) / 2;
    
    return xi_param
end


function nrbmak(coefs,knots)
    np = size(coefs);
    dim = np[1];
    number = np[2];
    if (dim < 4)
        # nurbs.coefs = repmat([0.0 0.0 1.0]',[1 np[2]]);
        # nurbs.coefs[1:dim,:] = coefs;
        coefs=[coefs;ones(1,np[2])]
    end
    order = length(knots)-np[2];
    knots = sort(vec(knots));
    knots = (knots-knots[1])/(knots[end]-knots[1]);
    
    nurbs=Curve(knots,coefs,number,order)
end

function nrbeval(nurbs,tt,basis=0)
    # Function Name:
    #
    #   nrbeval - Evaluate a NURBS at parameteric points
    #
    # Calling Sequence:
    #
    #   p = nrbeval(crv,ut)
    #   [p,w] = nrbeval(srf,{ut,vt})
    #
    # Parameters:
    #
    #   crv		: NURBS curve, see nrbmak.
    #
    #   srf		: NURBS surface, see nrbmak.
    #
    #   ut		: Parametric evaluation points along U direction.
    #
    #   vt		: Parametric evaluation points along V direction.
    #
    #   p		: Evaluated points on the NURBS curve or surface as cartesian
    # 		coordinates (x,y,z). If w is included on the lhs argument list
    # 		the points are returned as homogeneous coordinates (wx,wy,wz).
    #
    #   w		: Weights of the homogeneous coordinates of the evaluated
    # 		points. Note inclusion of this argument changes the type
    # 		of coordinates returned in p (see above).
    
    # NURBS structure represents a curve
    #  tt represent a vector of parametric points in the u direction
    
    val = bspeval(nurbs.order-1,nurbs.coefs,nurbs.knots,tt);
    w = val[4,:];
    p = val[1:3,:];
    p = p./repmat(w',3,1);
    p
end

function  bspeval(d,c,k,u)
    #
    # Function Name:
    #
    #   bspeval - Evaluate a univariate B-Spline.
    #
    # Calling Sequence:
    #
    #   p = bspeval(d,c,k,u)
    #
    # Parameters:
    #
    #   d	: Degree of the B-Spline.
    #
    #   c	: Control Points, matrix of size (dim,nc).
    #
    #   k	: Knot sequence, row vector of size nk.
    #
    #   u	: Parametric evaluation points, row vector of size nu.
    #
    #   p	: Evaluated points, matrix of size (dim,nu)
    #
    # Description:
    #
    #   Evaluate a univariate B-Spline. This function provides an interface to
    #   a toolbox 'C' routine.
    nu = length(u)
    mc,nc = size(c)
    #   int bspeval(int d, double *c, int mc, int nc, double *k, int nk, double *u,int nu, double *p){
    #   int ierr = 0;
    #   int i, s, tmp1, row, col;
    #   double tmp2;
    #
    #   // Construct the control points
    #   double **ctrl = vec2mat(c,mc,nc);
    #
    #   // Contruct the evaluated points
    p = zeros(mc,nu);                               #   double **pnt = vec2mat(p,mc,nu);
    #
    #   // space for the basis functions
    N = zeros(d+1,1);                               #   double *N = (double*) mxMalloc((d+1)*sizeof(double));
    #
    #   // for each parametric point i
    for col=1:nu                                    #   for (col = 0; col < nu; col++) {
        s = findspan(nc-1, d, u[col], k);           #     s = findspan(nc-1, d, u[col], k);
        N = basisfun(s,u[col],d,k);                 #     basisfun(s, u[col], d, k, N);
        #
        tmp1 = s - d + 1;                           #     tmp1 = s - d;
        for row=1:mc                                #     for (row = 0; row < mc; row++)  {
            tmp2 = 0;                               #       tmp2 = 0.0;
            for i=0:d                               #       for (i = 0; i <= d; i++)
                tmp2 = tmp2 + N[i+1]*c[row,tmp1+i];  # 	tmp2 += N[i] * ctrl[tmp1+i][row];
            end                                     #
            p[row,col] = tmp2;                      #       pnt[col][row] = tmp2;
        end                                         #     }
    end                                             #   }
    #
    #   mxFree(N);
    #   freevec2mat(pnt);
    #   freevec2mat(ctrl);
    #
    #   return ierr;
    p                                                #   }
end
function  findspan(n,p,u,U)
    # FINDSPAN  Find the span of a B-Spline knot vector at a parametric point
    # -------------------------------------------------------------------------
    # ADAPTATION of FINDSPAN from C
    # -------------------------------------------------------------------------
    #
    # Calling Sequence:
    #
    #   s = findspan(n,p,u,U)
    #
    #  INPUT:
    #
    #    n - number of control points - 1
    #    p - spline degree
    #    u - parametric point
    #    U - knot sequence
    #
    #  RETURN:
    #
    #    s - knot span
    #
    #  Algorithm A2.1 from 'The NURBS BOOK' pg68
    
    # int findspan(int n, int p, double u, double *U) {
    
    #   int low, high, mid;
    #   // special case
    if u==U[n+2]
        s=n
        return s
    end               #   if (u == U[n+1]) return(n);
    #
    #   // do binary search
    low = p;                                        #   low = p;
    high = n + 1;                                   #   high = n + 1;
    mid = floor(Int,(low + high) / 2);                  #   mid = (low + high) / 2;
    while (u < U[mid+1] || u >= U[mid+2])           #   while (u < U[mid] || u >= U[mid+1])  {
        if (u < U[mid+1])                           #     if (u < U[mid])
            high = mid;                             #       high = mid;
        else                                        #     else
            low = mid;                              #       low = mid;
        end
        mid = floor(Int,(low + high) / 2);              #     mid = (low + high) / 2;
    end                                             #   }
    #
    s = mid;                                        #   return(mid);
    s                                           #   }
end

function basisfun(i,u,p,U)
    # BASISFUN  Basis function for B-Spline
    # -------------------------------------------------------------------------
    # ADAPTATION of BASISFUN from C Routine
    # -------------------------------------------------------------------------
    #
    # Calling Sequence:
    #
    #   N = basisfun(i,u,p,U)
    #
    #    INPUT:
    #
    #      i - knot span  ( from FindSpan() )
    #      u - parametric point
    #      p - spline degree
    #      U - knot sequence
    #
    #    OUTPUT:
    #
    #      N - Basis functions vector[p+1]
    #
    #    Algorithm A2.2 from 'The NURBS BOOK' pg70.
    
    #   void basisfun(int i, double u, int p, double *U, double *N) {
    #   int j,r;
    #   double saved, temp;
    i = i + 1;
    #   // work space
    left = zeros(p+1,1);                              #   double *left  = (double*) mxMalloc((p+1)*sizeof(double));
    right = zeros(p+1,1);                             #   double *right = (double*) mxMalloc((p+1)*sizeof(double));
    N = zeros(p+1,1);
    N[1] = 1;                                         #   N[0] = 1.0;
    for j=1:p                                         #   for (j = 1; j <= p; j++) {
        left[j+1] = u - U[i+1-j];                     #   left[j]  = u - U[i+1-j];
        right[j+1] = U[i+j] - u;                      #   right[j] = U[i+j] - u;
        saved = 0;                                    #   saved = 0.0;
        
        for r=0:j-1                                   #   for (r = 0; r < j; r++) {
            temp = N[r+1]/(right[r+2] + left[j-r+1]); #   temp = N[r] / (right[r+1] + left[j-r]);
            N[r+1] = saved + right[r+2]*temp;         #   N[r] = saved + right[r+1] * temp;
            saved = left[j-r+1]*temp;                 #   saved = left[j-r] * temp;
        end                                           #   }
        
        N[j+1] = saved;                               #   N[j] = saved;
    end                                               #   }
    
    #   mxFree(left);
    #   mxFree(right);
    #   }
    N
end

function  bspline_deriv(order, knots, ctrl)
    # Knots and control points associated with the derivative of B-spline curve.
    #
    # Input arguments:
    # order:
    #    B-spline order (2 for linear, 3 for quadratic, etc.)
    # knots:
    #    knot vector
    # ctrl:
    #    control points, typically 2-by-m, 3-by-m, or 4-by-m (for weights)
    #
    # Output arguments:
    # dctrl:
    #    control points of the derivative of the input B-spline curve
    # dknots:
    #    the new knot vector associated with the derivative B-spline curve
    
    
    p = order - 1;
    tmp = size(ctrl);
    n = tmp[2]-1;
    dim = tmp[1];
    
    # derivative knots
    dknots = knots[2:end-1];
    
    # derivative control points
    dctrl = zeros(dim,n);
    for i = 1 : n
        dctrl[:,i] = (p / (knots[i+p+1] - knots[i+1])) * (ctrl[:,i+1] - ctrl[:,i]);
    end
    
    dctrl,dknots
end

function  nrbderiv(nurbs)
    #
    # Function Name:
    #
    #   nrbderiv - Construct the first derivative representation of a
    #              NURBS curve or surface.
    
    # NURBS structure represents a curve
    dcoefs,dknots = bspline_deriv(nurbs.order,nurbs.knots,nurbs.coefs);
    dnurbs = nrbmak(dcoefs, dknots);
end

function nrbdeval(nurbs, dnurbs, tt)
    val = bspeval(nurbs.order-1,nurbs.coefs,nurbs.knots,tt);
    cw = val[4,:];
    cp = val[1:3,:];
    # NURBS is a curve
    temp = repmat(cw,3,1);
    pnt = cp./temp;
    
    
    # first derivative
    #[cup,cuw] = nrbeval(dnurbs,tt);
    val = bspeval(dnurbs.order-1,dnurbs.coefs,dnurbs.knots,tt)
    cuw = val[4,:];
    cup = val[1:3,:];
    temp1 = repmat(cuw,3,1);
    jac = (cup-temp1.*pnt)./temp;
    pnt,jac
end

function nrbbasisfun(nrb,points)
    n    = size(nrb.coefs, 2) -1;
    p    = nrb.order -1;
    U    = nrb.knots;
    w    = nrb.coefs[4,:];
    npoints=length(points)
    B=zeros(nrb.order,npoints)
    nbfu=zeros(Int,nrb.order,npoints)
    for i=1:npoints
        
        
        spu  =  findspan(n, p, points[i], U);
        nbfu[:,i] =  (spu-p+1)+collect(0:p)
        
        N     = w[nbfu[:,i]] .* basisfun(spu, points[i], p, U);
        B[:,i]     = N./sum(N);
    end
    B, nbfu
end

function nrbplot(crv,divs=100)
    p=nrbeval(crv,linspace(0,1,divs))
    f=plot(p[1,:],p[2,:],c=:black,linewidth = 3)
    f=scatter(crv.coefs[1,:]./crv.coefs[3,:],crv.coefs[2,:]./crv.coefs[3,:],c=:black)
end

function bspdegelev(d,c,k,t)
    
    # BSPDEGELEV:  Degree elevate a univariate B-Spline.
    #
    # Calling Sequence:
    #
    #   [ic,ik] = bspdegelev(d,c,k,t)
    #
    # Parameters:
    #
    #   d - Degree of the B-Spline.
    #   c - Control points, matrix of size (dim,nc).
    #   k - Knot sequence, row vector of size nk.
    #   t - Raise the B-Spline degree t times.
    #
    #   ic - Control points of the new B-Spline.
    #   ik - Knot vector of the new B-Spline.
    #
    
    # O número de nós (knots), nknot, o número de pontos de controle, nc*t+2, e a ordem da
    # curva, d+t , estão relacionados por:
    #                           nknot = (nc*t+2) + (d+t) + 1
    mc,nc = size(c);
    # nc = number of control points
    # nknot = number of knots
    nknot=length(k)
    #println(nc*(t+1)+1)
    #println(nknot*(t+1)+1)
    ic = zeros(mc,nc);                                  #   double **ictrl = vec2mat(ic, mc, nc*(t+1));
    n = nc - 1;                                               #   n = nc - 1;
    #
    bezalfs =  zeros(d+1,d+t+1);                              #   bezalfs = matrix(d+1,d+t+1);
    bpts = zeros(mc,d+1);                                     #   bpts = matrix(mc,d+1);
    ebpts = zeros(mc,d+t+1);                                  #   ebpts = matrix(mc,d+t+1);
    Nextbpts = zeros(mc,d+1);                                 #   Nextbpts = matrix(mc,d+1);
    alfs = zeros(d,1);                                        #   alfs = (double *) mxMalloc(d*sizeof(double));
    #
    m = n + d + 1;                                            #   m = n + d + 1;
    ph = d + t;                                               #   ph = d + t;
    ph2 = Int32(floor(ph / 2));                                      #   ph2 = ph / 2;
    ik=zeros(ph+1)  # Duvida !!!!!!!!!!!!
    
    #
    #   // compute bezier degree elevation coefficeients
    bezalfs[1,1] = 1;                                         #   bezalfs[0][0] = bezalfs[ph][d] = 1.0;
    bezalfs[d+1,ph+1] = 1;                                    #
    
    for i=1:1:ph2
        inv = 1/bincoeff(ph,i)
        #     inv = 1.0 / bincoeff(ph,i);
        mpi = minimum([d,i]);                                        #     mpi = min(d,i);
        #
        for j=maximum([0,i-t]):mpi                                   #     for (j = max(0,i-t); j <= mpi; j++)
            bezalfs[j+1,i+1] = inv*bincoeff(d,j)*bincoeff(t,i-j)  #       bezalfs[i][j] = inv * bincoeff(d,j) * bincoeff(t,i-j);
        end
    end                                                       #   }
    #
    for i=ph2+1:ph-1                                          #   for (i = ph2+1; i <= ph-1; i++) {
        mpi = minimum([d,i]);                                        #     mpi = min(d, i);
        for j=maximum([0,i-t]):mpi                                   #     for (j = max(0,i-t); j <= mpi; j++)
            bezalfs[j+1,i+1] = bezalfs[d-j+1,ph-i+1];          #       bezalfs[i][j] = bezalfs[ph-i][d-j];
        end
    end                                                       #   }
    #
    mh = ph;                                                  #   mh = ph;
    kind = ph+1;                                              #   kind = ph+1;
    r = -1;                                                   #   r = -1;
    a = d;                                                    #   a = d;
    b = d+1;                                                  #   b = d+1;
    cind = 1;                                                 #   cind = 1;
    ua = k[1];                                                #   ua = k[0];
    #
    for ii=0:mc-1                                             #   for (ii = 0; ii < mc; ii++)
        ic[ii+1,1] = c[ii+1,1]                                #     ictrl[0][ii] = ctrl[0][ii];
    end                                                       #
    for i=0:ph                                                #   for (i = 0; i <= ph; i++)
        ik[i+1] = ua;                                          #     ik[i] = ua;
    end                                                       #
    #   // initialise first bezier seg
    for i=0:d                                                 #   for (i = 0; i <= d; i++)
        for ii=0:mc-1                                          #     for (ii = 0; ii < mc; ii++)
            bpts[ii+1,i+1] = c[ii+1,i+1];                       #       bpts[i][ii] = ctrl[i][ii];
        end
    end                                                       #
    #   // big loop thru knot vector
    while b < m                                               #   while (b < m)  {
        i = b;                                                 #     i = b;
        while b < m && k[b+1] == k[b+2]                        #     while (b < m && k[b] == k[b+1])
            b = b + 1;                                          #       b++;
        end                                                    #
        mul = b - i + 1;                                       #     mul = b - i + 1;
        mh = mh + mul + t;                                     #     mh += mul + t;
        ub = k[b+1];                                           #     ub = k[b];
        oldr = r;                                              #     oldr = r;
        r = d - mul;                                           #     r = d - mul;
        #
        #     // insert knot u(b) r times
        if oldr > 0                                            #     if (oldr > 0)
            lbz = floor((oldr+2)/2);                            #       lbz = (oldr+2) / 2;
        else                                                   #     else
            lbz = 1;                                            #       lbz = 1;
        end                                                    #
        
        if r > 0                                               #     if (r > 0)
            rbz = ph - floor((r+1)/2);                          #       rbz = ph - (r+1)/2;
        else                                                   #     else
            rbz = ph;                                           #       rbz = ph;
        end                                                    #
        
        if r > 0                                               #     if (r > 0) {
            #       // insert knot to get bezier segment
            numer = ub - ua;                                    #       numer = ub - ua;
            for q=d:-1:mul+1                                    #       for (q = d; q > mul; q--)
                alfs[q-mul] = numer / (k[a+q+1]-ua);             #         alfs[q-mul-1] = numer / (k[a+q]-ua);
            end
            
            for j=1:r                                           #       for (j = 1; j <= r; j++)  {
                save = r - j;                                    #         save = r - j;
                s = mul + j;                                     #         s = mul + j;
                #
                for q=d:-1:s                                     #         for (q = d; q >= s; q--)
                    for ii=0:mc-1                                 #           for (ii = 0; ii < mc; ii++)
                        tmp1 = alfs[q-s+1]*bpts[ii+1,q+1]
                        tmp2 = (1-alfs[q-s+1])*bpts[ii+1,q]
                        bpts[ii+1,q+1] = tmp1 + tmp2;              #             bpts[q][ii] = alfs[q-s]*bpts[q][ii]+(1.0-alfs[q-s])*bpts[q-1][ii];
                    end
                end                                              #
                
                for ii=0:mc-1                                    #         for (ii = 0; ii < mc; ii++)
                    Nextbpts[ii+1,save+1] = bpts[ii+1,d+1]       #           Nextbpts[save][ii] = bpts[d][ii];
                end
            end                                                 #       }
        end                                                    #     }
        #     // end of insert knot
        #
        #     // degree elevate bezier
        for i=lbz:ph                                           #     for (i = lbz; i <= ph; i++)  {
            for ii=0:mc-1                                       #       for (ii = 0; ii < mc; ii++)
                ebpts[ii+1,i+1] = 0;                             #         ebpts[i][ii] = 0.0;
            end
            mpi = minimum([d, i]);                                    #       mpi = min(d, i);
            for j=maximum([0,i-t]):mpi                                #       for (j = max(0,i-t); j <= mpi; j++)
                for ii=0:mc-1                                    #         for (ii = 0; ii < mc; ii++)
                    tmp1 = ebpts[ii+1,i+1]
                    tmp2 = bezalfs[j+1,i+1]*bpts[ii+1,j+1]
                    ebpts[ii+1,i+1] = tmp1 + tmp2;                #           ebpts[i][ii] = ebpts[i][ii] + bezalfs[i][j]*bpts[j][ii];
                end
            end
        end                                                    #     }
        #     // end of degree elevating bezier
        #
        if oldr > 1                                            #     if (oldr > 1)  {
            #       // must remove knot u=k[a] oldr times
            first = kind - 2;                                                    #       first = kind - 2;
            last = kind;                                        #       last = kind;
            den = ub - ua;                                      #       den = ub - ua;
            bet = floor((ub-ik[kind]) / den);                   #       bet = (ub-ik[kind-1]) / den;
            #
            #       // knot removal loop
            for tr=1:oldr-1                                     #       for (tr = 1; tr < oldr; tr++)  {
                i = first;                                       #         i = first;
                j = last;                                        #         j = last;
                kj = j - kind + 1;                               #         kj = j - kind + 1;
                while j-i > tr                                   #         while (j - i > tr)  {
                    #           // loop and compute the new control points
                    #           // for one removal step
                    if i < cind                                   #           if (i < cind)  {
                        alf = (ub-ik[i+1])/(ua-ik[i+1]);           #             alf = (ub-ik[i])/(ua-ik[i]);
                        nic=size(ic,2)                          #             for (ii = 0; ii < mc; ii++)
                        tmp1 = alf*ic[:,i+1]
                        tmp2 = (1-alf)*ic[:,i]
			if(nic>=i+1)
                            ic[:,i+1] = tmp1 + tmp2;
                        else
                            ic=[ic (tmp1+tmp2)]					   #               ictrl[i][ii] = alf * ictrl[i][ii] + (1.0-alf) * ictrl[i-1][ii];
                        end
                    end                                           #           }
                    if j >= lbz                                   #           if (j >= lbz)  {
                        if j-tr <= kind-ph+oldr                    #             if (j-tr <= kind-ph+oldr) {
                            gam = (ub-ik[j-tr+1]) / den;            #               gam = (ub-ik[j-tr]) / den;
                            for ii=0:mc-1                           #               for (ii = 0; ii < mc; ii++)
                                tmp1 = gam*ebpts[ii+1,kj+1]
                                tmp2 = (1-gam)*ebpts[ii+1,kj+2]
                                ebpts[ii+1,kj+1] = tmp1 + tmp2;      #                 ebpts[kj][ii] = gam*ebpts[kj][ii] + (1.0-gam)*ebpts[kj+1][ii];
                            end                                     #             }
                        else                                       #             else  {
                            for ii=0:mc-1                           #               for (ii = 0; ii < mc; ii++)
                                tmp1 = bet*ebpts[ii+1,kj+1];
                                tmp2 = (1-bet)*ebpts[ii+1,kj+2];
                                ebpts[ii+1,kj+1] = tmp1 + tmp2;      #                 ebpts[kj][ii] = bet*ebpts[kj][ii] + (1.0-bet)*ebpts[kj+1][ii];
                            end
                        end                                        #             }
                    end                                           #           }
                    i = i + 1;                                    #           i++;
                    j = j - 1;                                    #           j--;
                    kj = kj - 1;                                  #           kj--;
                end                                              #         }
                #
                first = first - 1;                               #         first--;
                last = last + 1;                                 #         last++;
            end                                                 #       }
        end                                                    #     }
#     // end of removing knot n=k[a]
#
#     // load the knot ua
if a != d                                              #     if (a != d)
    for i=0:ph-oldr-1                                   #       for (i = 0; i < ph-oldr; i++)  {
	nik=length(ik)
	if(nik>=kind+1)
            ik[kind+1] = ua;                                 #         ik[kind] = ua;
        else
	    append!(ik,ua)
	end
	kind = kind + 1;                                 #         kind++;
    end
end                                                    #       }
#
#     // load ctrl pts int
for j=lbz:rbz                                       #     for (j = lbz; j <= rbz; j++)  {
    nic=size(ic,2)
    if(nic>=cind+1)
	ic[:,cind+1] = ebpts[:,j+1]            #         ictrl[cind][ii] = ebpts[j][ii];
    else
        ic=[ic ebpts[:,j+1]]
    end
    cind = cind + 1;                                 #       cind++;
end                                                 #     }
#
if b < m                                            #     if (b < m)  {
    #       // setup for next pass thru loop
    for j=0:r-1                                      #       for (j = 0; j < r; j++)
        for ii=0:mc-1                                 #         for (ii = 0; ii < mc; ii++)
            bpts[ii+1,j+1] = Nextbpts[ii+1,j+1]       #           bpts[j][ii] = Nextbpts[j][ii];
        end
    end
    for j=r:d                                        #       for (j = r; j <= d; j++)
        for ii=0:mc-1                                 #         for (ii = 0; ii < mc; ii++)
            bpts[ii+1,j+1] = c[ii+1,b-d+j+1]          #           bpts[j][ii] = ctrl[b-d+j][ii];
        end
    end
    a = b;                                           #       a = b;
    b = b+1;                                         #       b++;
    ua = ub;                                         #       ua = ub;
    #     }
else                                                #     else
    #       // end knot
    for i=0:ph
        nik=length(ik)
	if(nik>=kind+i+1)
	    ik[kind+i+1] = ub;                            #         ik[kind+i] = ub;
        else
	    append!(ik,ub)
	end
    end
end
end
return ic,ik
end



function bincoeff(n,k)
    #  Computes the binomial coefficient.
    #
    #      ( n )      n!
    #      (   ) = --------
    #      ( k )   k!(n-k)!
    #
    #  b = bincoeff(n,k)
    #
    #  Algorithm from 'Numerical Recipes in C, 2nd Edition' pg215.
    
    # double bincoeff(int n, int k)
    # {
    b = floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));      #   return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
    # }
    return b
end

function factln(n)
    # computes ln(n!)
    if n <= 1
        f = 0;
    else
        f = log(factorial(n)) #log(factorial(n));
    end
    return f
end



function  bspkntins(d,c,k,u)
    
    # BSPKNTINS:  Insert knots into a B-Spline
    #
    # Calling Sequence:
    #
    #   [ic,ik] = bspkntins(d,c,k,u)
    #
    #  INPUT:
    #
    #    d - spline degree             integer
    #    c - control points            double  matrix(mc,nc)
    #    k - knot sequence             double  vector(nk)
    #    u - new knots                 double  vector(nu)
    #
    #  OUTPUT:
    #
    #    ic - new control points double  matrix(mc,nc+nu)
    #    ik - new knot sequence  double  vector(nk+nu)
    #
    #  Modified version of Algorithm A5.4 from 'The NURBS BOOK' pg164.
    #
    
    mc,nc = size(c);
    nu = length(u);
    nk = length(k);
    #
    # int bspkntins(int d, double *c, int mc, int nc, double *k, int nk,
    #               double *u, int nu, double *ic, double *ik)
    # {
    #   int ierr = 0;
    #   int a, b, r, l, i, j, m, n, s, q, ind;
    #   double alfa;
    #
    #   double **ctrl  = vec2mat(c, mc, nc);
    ic = zeros(mc,nc+nu);                                #   double **ictrl = vec2mat(ic, mc, nc+nu);
    ik = zeros(nk+nu);
    #
    n = size(c,2) - 1;                                   #   n = nc - 1;
    r = length(u) - 1;                                   #   r = nu - 1;
    #
    m = n + d + 1;                                       #   m = n + d + 1;
    a = findspan(n, d, u[1], k);                         #   a = findspan(n, d, u[0], k);
    b = findspan(n, d, u[r+1], k);                       #   b = findspan(n, d, u[r], k);
    b = b+1;                                             #   ++b;
    #
    for q=0:mc-1                                         #   for (q = 0; q < mc; q++)  {
        for j=0:a-d
            ic[q+1,j+1] = c[q+1,j+1]
        end        #     for (j = 0; j <= a-d; j++) ictrl[j][q] = ctrl[j][q];
        for j=b-1:n
            ic[q+1,j+r+2] = c[q+1,j+1]
        end      #     for (j = b-1; j <= n; j++) ictrl[j+r+1][q] = ctrl[j][q];
    end                                                  #   }
    
    for j=0:a
        ik[j+1] = k[j+1]
    end                     #   for (j = 0; j <= a; j++)   ik[j] = k[j];
    for j=b+d:m
        ik[j+r+2] = k[j+1]
    end                 #   for (j = b+d; j <= m; j++) ik[j+r+1] = k[j];
    #
    i = b + d - 1;                                       #   i = b + d - 1;
    s = b + d + r;                                       #   s = b + d + r;
    
    for j=r:-1:0                                         #   for (j = r; j >= 0; j--) {
        while u[j+1] <= k[i+1] && i > a                   #     while (u[j] <= k[i] && i > a) {
            for q=0:mc-1                                  #       for (q = 0; q < mc; q++)
                ic[q+1,s-d] = c[q+1,i-d]                 #         ictrl[s-d-1][q] = ctrl[i-d-1][q];
            end
            ik[s+1] = k[i+1]                             #       ik[s] = k[i];
            s = s - 1;                                    #       --s;
            i = i - 1;                                    #       --i;
        end                                               #     }
        
        for q=0:mc-1                                      #     for (q = 0; q < mc; q++)
            ic[q+1,s-d] = ic[q+1,s-d+1]                  #       ictrl[s-d-1][q] = ictrl[s-d][q];
        end
        
        for l=1:d                                         #     for (l = 1; l <= d; l++)  {
            ind = s - d + l;                              #       ind = s - d + l;
            alfa = ik[s+l+1] - u[j+1]                    #       alfa = ik[s+l] - u[j];
            if abs(alfa) == 0                             #       if (fabs(alfa) == 0.0)
                for q=0:mc-1                              #         for (q = 0; q < mc; q++)
                    ic[q+1,ind] = ic[q+1,ind+1]          #           ictrl[ind-1][q] = ictrl[ind][q];
                end
            else                                          #       else  {
                alfa = alfa/(ik[s+l+1] - k[i-d+l+1]);     #         alfa /= (ik[s+l] - k[i-d+l]);
                for q=0:mc-1                              #         for (q = 0; q < mc; q++)
                    tmp = (1-alfa)*ic[q+1,ind+1];
                    ic[q+1,ind] = alfa*ic[q+1,ind] + tmp; #           ictrl[ind-1][q] = alfa*ictrl[ind-1][q]+(1.0-alfa)*ictrl[ind][q];
                end
            end                                           #       }
        end                                               #     }
        #
        ik[s+1] = u[j+1]                                 #     ik[s] = u[j];
        s = s - 1;                                        #     --s;
    end
    return ic,ik
end

function nrbline(p1,p2)
    #
    # Function Name:
    #
    #   nrbline - Construct a straight line.
    #
    # Calling Sequence:
    #
    #   crv = nrbline()
    #   crv = nrbline(p1,p2)
    #
    # Parameters:
    #
    # p1		: 2D or 3D cartesian coordinate of the start point.
    #
    # p2            : 2D or 3D cartesian coordinate of the end point.
    #
    # crv		: NURBS curve for a straight line.
    #
    # Description:
    #
    #   Constructs NURBS data structure for a straight line. If no rhs
    #   coordinates are included the function returns a unit straight
    #   line along the x-axis.
    
    #  D.M. Spink
    #  Copyright (c) 2000.
    
    coefs = [zeros(3,2); ones(1,2)];
    
    coefs[1:length(p1),1] = p1[:];
    coefs[1:length(p2),2] = p2[:];
    
    curve = nrbmak(coefs, [0 0 1 1]);
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

    eta = A*gamm.^3 + B*gamm.^2 + C*gamm + D;
    Jt = 3*A*gamm.^2 + 2*B*gamm + C;
    return  eta,Jt
end
