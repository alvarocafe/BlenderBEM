## Test cases for wave problems
# This script tests the BEM models against known problems
#for which analytical solutions are readily obtained.
#include("../BEM_base.jl")
#################### Test case n ####################
### Problem description
### Analytical solution
### BEM model
### return error
 using Plots
#gr()
#pyplot()
xs = 0.01:0.01:1
# p1 = contour(xs,xs,f,fill=true) 
#################### Test case 1-4 ####################
### Acoustic tube
# Consider a square acoustic domain in which the speed of
#sound is approximately 344 [m/s] inside the domain and
#there are two opposite walls which are rigid and two
#open-ended, such that one of them is excited at a frequency
# ω [rads], so that the wave number is k = ω/c, with
#relative amplitude P = 1.
### Analytical solution
phi_closed(k,x) = cos.(k.*x./L)
q_closed(k,x) = -k.*sin(k.*x)
kclosed(L=1,n=1) = n*π./(L)

function nurbsclosed2D(ne=10,L=1,k=kclosed())
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    POINTS = [1 0 0
	      2 L 0
	      3 L L
	      4 0 L];
    SEGMENTS = [1 1 2 0
	        2 2 3 0
	        3 3 4 0
	        4 4 1 0];
    #PONTOS_dom = [1 L/2 L/2]
    n_pint = 100; # Number of domain points
    PONTOS_dom = zeros(n_pint,3);
    delta = 0.01; # distance from both ends 
    passo = (L-2*delta)/(n_pint-1);
    for i = 1:n_pint
        PONTOS_dom[i,:] = [i delta+(i-1)*passo L/2];
    end
    fc = [0 0 0]
    BCSeg = [1 1 0
	     2 1 0
	     3 1 0
	     4 0 1];
    info = [POINTS,SEGMENTS,BCSeg,k];
    t = @elapsed phi,q,phi_dom = wavenurbs2D.solveH(info, PONTOS_dom, fc, k)  
    return norm(phi_dom.^2 .- phi_closed(k,PONTOS_dom[1,2]).^2)./size(PONTOS_dom,1)
end

### BEM model
function constclosed2D(ne=10,L=1,k=kclosed())
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    POINTS = [1 0 0
	      2 L 0
	      3 L L
	      4 0 L];
    SEGMENTS = [1 1 2 0
	        2 2 3 0
	        3 3 4 0
	        4 4 1 0];
    # Each segment will be meshed by ne elements
    #ne = 100;
    MESH = [1 ne
	    2 ne
	    3 ne
	    4 ne];
    # PONTOS_dom = [1 L/2 L/2]
    n_pint = 100; # Number of domain points
    PONTOS_dom = zeros(n_pint,3);
    delta = 0.01; # distance from both ends 
    passo = (L-2*delta)/(n_pint-1);
    for i = 1:n_pint
        PONTOS_dom[i,:] = [i delta+(i-1)*passo L/2];
    end
    fc = [0 0 0]
    BCSeg = [1 1 0
	     2 1 0
	     3 1 0
	     4 1 1];
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    t = @elapsed phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    # Conventional method with approximated influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    #    info = [NOS_GEO,NOS,ELEM,CDC]
    #    tH = @elapsed phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)    
    return norm(phi_dom.^2 .- phi_closed(k,PONTOS_dom[1,2]).^2)./size(PONTOS_dom,1)
end

### Analytical solution
phi_cup(k,x) = sin.(k.*x)
q_cup(k,x) = -k.*cos(k.*x)
function cup2D(ne=10,L=1,k=1)
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    POINTS = [1 0 0
	      2 L 0
	      3 L L
	      4 0 L];
    SEGMENTS = [1 1 2 0
	        2 2 3 0
	        3 3 4 0
	        4 4 1 0];
    # Each segment will be meshed by ne elements
    #ne = 100;
    MESH = [1 ne
	    2 ne
	    3 ne
	    4 ne];
    PONTOS_dom = [1 L/2 L/2]
    fc = [0 0 0]
    BCSeg = [1 1 0
	     2 0 1
	     3 1 0
	     4 1 0];
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    t = @elapsed phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    # Conventional method with approximated influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    #    info = [NOS_GEO,NOS,ELEM,CDC]
    #    tH = @elapsed phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)    
    return norm(phi_dom-phi_cup(k,PONTOS_dom[1,2]))
end

### Analytical solution
phi_open(k,x) = sin.(k.*x)
q_open(k,x) = -k.*cos(k.*x)
function open2D(ne=10,L=1,k=1)
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    POINTS = [1 0 0
	      2 L 0
	      3 L L
	      4 0 L];
    SEGMENTS = [1 1 2 0
	        2 2 3 0
	        3 3 4 0
	        4 4 1 0];
    # Each segment will be meshed by ne elements
    #ne = 100;
    MESH = [1 ne
	    2 ne
	    3 ne
	    4 ne];
    PONTOS_dom = [1 L/2 L/2]
    fc = [0 0 0]
    BCSeg = [1 1 0
	     2 0 0
	     3 1 0
	     4 0 1];
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    t = @elapsed phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    # Conventional method with approximated influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    #    info = [NOS_GEO,NOS,ELEM,CDC]
    #    tH = @elapsed phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)    
    return norm(phi_dom-phi_open(k,PONTOS_dom[1,2]))
end

### Analytical solution
phi_helmholtz_bottle(k,x) = sin.(k.*x)
q_helmholtz_bottle(k,x) = -k.*cos(k.*x)
function helmholtz_bottle(ne=10,L=1,k=1)
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    POINTS = [1 0 0
	      2 L 0
	      3 L L
	      4 0 L];
    SEGMENTS = [1 1 2 0
	        2 2 3 0
	        3 3 4 0
	        4 4 1 0];
    # Each segment will be meshed by ne elements
    #ne = 100;
    MESH = [1 ne
	    2 ne
	    3 ne
	    4 ne];
    PONTOS_dom = [1 L/2 L/2]
    fc = [0 0 0]
    BCSeg = [1 1 0
	     2 0 0
	     3 1 0
	     4 0 1];
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    t = @elapsed phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    # Conventional method with approximated influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    #    info = [NOS_GEO,NOS,ELEM,CDC]
    #    tH = @elapsed phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)    
    return norm(phi_dom-phi_helmholtz(k,PONTOS_dom[1,2]))
end


#################### Test case 5 ####################
### Vibrating cylinder
# An acoustically rigid cylinder is let to vibrate in an
#infinite acoustic domain at frequency ω so that the wavenumber is k = ω/c.
### Analytical solution
#The analytical solution for the vibrating cylinder obtained
#by applying variable separation in cylindrical coordinates
#of the Helmholtz equation is
phi_cylinder(k,r,x) = (1/k).*(SpecialFunctions.besselh(0,2,k.*x)./SpecialFunctions.besselh(1,2,k.*r));
#at a distance x from the cylinder of radius r.
### BEM models
function nurbscylinder(r=0.5,c=[0 0],k=1)
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    POINTS =[1  c[1,1]-r	c[1,2]
       	     2	c[1,1]		c[1,2]+r
       	     3	c[1,1]+r	c[1,2]
       	     4	c[1,1]		c[1,2]-r];
    SEGMENTS = [1 1 2 -r
                2 2 3 -r
                3 3 4 -r
                4 4 1 -r];
    BCSeg = [1 1 1 0
             2 1 1 0
             3 1 1 0
             4 1 1 0];
    #PONTOS_dom = [1 L/2 L/2]
    n_pint = 100; # Number of domain points
    PONTOS_dom = zeros(n_pint,3);
    delta = 10*r; # distance from both ends 
    passo = delta/25
    for i = 1:n_pint
        PONTOS_dom[i,:] = [i delta+(i-1)*passo r];
    end
    fc = [0 0 0]
    BCSeg = [1 1 0
	     2 1 0
	     3 1 0
	     4 0 1];
    info = [POINTS,SEGMENTS,BCSeg,k];
    t = @elapsed phi,q,phi_dom = wavenurbs2D.solveH(info, PONTOS_dom, fc, k)  
    #return norm(phi_dom.^2 .- phi_cylinder(k,r,PONTOS_dom[:,2]).^2)./size(PONTOS_dom,1)
    return phi,q,phi_dom, PONTOS_dom
end


function const2Dcylinder(ne = 100,r=0.5,c=[0 0],k=1)
    t=0; ϵ=0; phi=0; q=0; phi_dom=0;  
    POINTS =[1  c[1,1]-r	c[1,2]
       	     2	c[1,1]		c[1,2]+r
       	     3	c[1,1]+r	c[1,2]
       	     4	c[1,1]		c[1,2]-r];
    SEGMENTS = [1 1 2 -r
                2 2 3 -r
                3 3 4 -r
                4 4 1 -r];
    BCSeg = [1 1 1 0
             2 1 1 0
             3 1 1 0
             4 1 1 0];
    MESH = [1 ne
            2 ne
            3 ne
            4 ne];
    x = 10*r
    PONTOS_dom = [1 0 x]
    fc = [0 0 0]
    # # Conventional method with full influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    t = @elapsed phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    b1 = 1:nnos # Array containing all the indexes for nodes and elements
    npg=6
    qsi,w = const2D.Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = const2D.cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    ϵ = abs.(sqrt.(((phi_dom .- phi_cylinder(k,r,x)).^2)./phi_cylinder(k,r,x).^2))
    return ϵ[1]
end

function Hconst2Dcylinder(ne = 100,r=0.5,c=[0 0],k=1)
    t=0; tH=0; tiso=0; ϵ=0; ϵH=0; ϵiso=0; phi=0; q=0; phi_dom=0;  phiH=0; qH=0; phi_domH=0;  phiiso=0; qiso=0; domiso=0;
    POINTS =[1  c[1,1]-r	c[1,2]
       	     2	c[1,1]		c[1,2]+r
       	     3	c[1,1]+r	c[1,2]
       	     4	c[1,1]		c[1,2]-r];
    SEGMENTS = [1 1 2 -r
                2 2 3 -r
                3 3 4 -r
                4 4 1 -r];
    BCSeg = [1 1 1 0
             2 1 1 0
             3 1 1 0
             4 1 1 0];
    MESH = [1 ne
            2 ne
            3 ne
            4 ne];
    x = 10*r
    PONTOS_dom = [1 0 x]
    fc = [0 0 0]
    # # Conventional method with full influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    # Conventional method with approximated influence matrices
    #    tH = @elapsed phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)
    Tree,block = const2D.cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
    tH = @elapsed Aaca1,baca =  const2D.Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    Aaca = const2D.montacheia(Aaca1,block,Tree,nnos)
    xaca = Aaca\baca
    phiH,qphiH,phidomH = const2D.(xaca,CDC)
    #println("A = ",A)
    #println("Aaca = ",Aaca)
    #ϵH = norm(Aaca-A)
    #println("ϵH = ",ϵH)
    ϵH = abs.(sqrt.(((phi_domH .- phi_cylinder(1,0.5,5)).^2)./phi_cylinder(1,0.5,5).^2))
    return ϵH[1]
end


function nurbs2Dcylinder(ne = 100,r=0.5,c=[0 0],k=1)
    t=0; tH=0; tiso=0; ϵ=0; ϵH=0; ϵiso=0; phi=0; q=0; phi_dom=0;  phiH=0; qH=0; phi_domH=0;  phiiso=0; qiso=0; domiso=0;
    POINTS =[1  c[1,1]-r	c[1,2]
       	     2	c[1,1]		c[1,2]+r
       	     3	c[1,1]+r	c[1,2]
       	     4	c[1,1]		c[1,2]-r];
    SEGMENTS = [1 1 2 -r
                2 2 3 -r
                3 3 4 -r
                4 4 1 -r];
    BCSeg = [1 1 1 0
             2 1 1 0
             3 1 1 0
             4 1 1 0];
    MESH = [1 ne
            2 ne
            3 ne
            4 ne];
    x = 10*r
    PONTOS_dom = [1 0 x]
    fc = [0 0 0]
    # Isogeometric BEM with full influence matrices
    collocCoord,nnos,crv,dcrv,E,CDC = nurbs2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    info = [collocCoord,nnos,crv,dcrv,E]
    tiso = @elapsed phiiso, qiso, phi_domiso = nurbs2D.solve(info,PONTOS_dom,[0],CDC,k)
    ϵ = abs.(sqrt.(((phi_domiso .- phi_cylinder(1,0.5,5)).^2)./phi_cylinder(1,0.5,5).^2))
    return ϵ[1]
end


function cylinder3D(ne = 100,r=1,L=10,c=[0 0],k=1)
    t=0; tH=0; tiso=0; ϵ=0; ϵH=0; ϵiso=0; phi=0; q=0; phi_dom=0;  phiH=0; qH=0; phi_domH=0;  phiiso=0; qiso=0; domiso=0;
    #3D constant elements    
    #NURBS
    PNTS = [1 0 0 0
            2 r 0 0
            3 r L 0
            4 0 L 0]
    SEGS = [1 1 2 0
            2 2 3 0
            3 3 4 0
            4 4 1 0]
    MESH = [1 10
            2 10
            3 10
            4 10]    
    BCSeg = [1 1 1
             2 1 1
             3 1 1
             4 1 1]
    crv = nurbs2D.format_dad(PNTS,SEGS,MESH,BCSeg)
    sup = Array{Surf}(1)
    sup[1] = nurbs2D.nrbrevolve(crv,[0.5, 0.5, 0.5],[1,0,0])
    return ϵ
end

#################### Test case 6 ####################
### Plane wave scattering on a rigid cylinder
# An acoustically rigid cylinder is let to vibrate in an
#infinite acoustic domain at frequency ω so that the wavenumber is k = ω/c.
# The cylinder also scatters a plane wave with direction d = (dx,dy,dz) and
#wavenumber k_inc, amplitude a_inc. 
### Analytical solution
#The analytical solution for the scattering cylinder is
besseljprime(nu,x) = -1/2*(SpecialFunctions.besselj(nu-1,x) - SpecialFunctions.besselj(nu+1,x))
besselyprime(nu,x) = -1/2*(SpecialFunctions.bessely(nu-1,x) - SpecialFunctions.bessely(nu+1,x))
besselhprime(nu,x) = besseljprime(nu,x) + complex(0,1)*besselyprime(nu,x)

function phi_cyl_scat(k,r,θ,a)
    phi = -(SpecialFunctions.besselj(1,k*a)/SpecialFunctions.besselh(1,k*a))*SpecialFunctions.besselh(0,k*r)
    for n in 1:1:100
        phi += (complex(0,1)^n)*(besseljprime(n,k*a)/besselhprime(n,k*a))*besselh(n,k*r)*cos(n*θ)
    end
    return phi
end
### BEM model
function cyl_scat(ne = 100,r=0.5,c=[0 0],k=1)
    t=0; tH=0; tiso=0; ϵ=0; ϵH=0; ϵiso=0; phi=0; q=0; phi_dom=0;  phiH=0; qH=0; phi_domH=0;  phiiso=0; qiso=0; domiso=0;
    POINTS =[1  c[1,1]-r	c[1,2]
       	     2	c[1,1]		c[1,2]+r
       	     3	c[1,1]+r	c[1,2]
       	     4	c[1,1]		c[1,2]-r];
    SEGMENTS = [1 1 2 -r
                2 2 3 -r
                3 3 4 -r
                4 4 1 -r];
    BCSeg = [1 1 1 0
             2 1 1 0
             3 1 1 0
             4 1 1 0];
    MESH = [1 ne
            2 ne
            3 ne
            4 ne];
    x = 10*r
    PONTOS_dom = [1 0 x]
    fc = [1 -1 0]
    # # Conventional method with full influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    #    t = @elapsed phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    b1 = 1:nnos # Array containing all the indexes for nodes and elements
    npg=6
    qsi,w = const2D.Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    t = @elapsed A,b = const2D.cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    ϵ = abs.(sqrt.(((phi_dom .- phi_cyl_scat(k,r,x)).^2)./phi_cyl_scat(k,r,x).^2))
    Tree,block = const2D.cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
    tH = @elapsed Aaca1,b =  const2D.Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    println(typeof(Aaca1))
    println(size(Aaca1[1]))
    println(size(Aaca1[2]))
    Aaca = const2D.montacheia(Aaca1,block,Tree,nnos)
    println("A = ",A)
    #println("Aaca = ",Aaca)
    ϵH = norm(Aaca-A)
    println("ϵH = ",ϵH)
    # Conventional method with approximated influence matrices
    #    tH = @elapsed phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)
    #    ϵH = abs.(sqrt.(((phi_domH .- phi_cylinder(1,0.5,5)).^2)./phi_cylinder(1,0.5,5).^2))
    # Isogeometric BEM with full influence matrices
    #collocCoord,nnos,crv,dcrv,E,CDC = nurbs2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    #info = [collocCoord,nnos,crv,dcrv,E]
    #tiso = @elapsed phiiso, qiso, phi_domiso = nurbs2D.solve(info,PONTOS_dom,[0],CDC,k)
    #ϵiso = abs.(sqrt.(((phi_domiso .- phi_cylinder(1,0.5,5)).^2)./phi_cylinder(1,0.5,5).^2))
    return ϵ
end

function cyl_scat3D(ne = 100,r=0.5,c=[0 0],k=1)
    t=0; tH=0; tiso=0; ϵ=0; ϵH=0; ϵiso=0; phi=0; q=0; phi_dom=0;  phiH=0; qH=0; phi_domH=0;  phiiso=0; qiso=0; domiso=0;
    return ϵH
end

#################### Test case 7 ####################
### Pulsating sphere
# An acoustically rigid sphere is let to vibrate in an
#infinite acoustic domain at frequency ω so that the wavenumber is k = ω/c.
### Analytical solution
#The analytical solution for the scattering cylinder is
phi_sphere(k,r,a,ρ,c) = a/r*ρ*c*(-complex(0,1)*k*a/(1+complex(0,1)*k*a))*exp(complex(0,1)*k*(r-a));
#at a distance x from the cylinder of radius r.
### BEM model

#################### Test case 8 ####################
### Plane wave scattering on a rigid sphere
# An acoustically rigid sphere is let to vibrate in an
#infinite acoustic domain at frequency ω so that the wavenumber is k = ω/c.
# The sphere also scatters a plane wave with direction d = (dx,dy,dz) and
#wavenumber k_inc, amplitude a_inc. 
### Analytical solution
#The analytical solution for the scattering cylinder is
function phi_sph_scat(k,r,θ,a)
    phi = complex(0,0)
    for n in 1:10
        phi += ((complex(0,1)^(n) * (2*n +1) * SpecialFunctions.bessel(n,1,k*a))/SpecialFunctions.hankel(n,1,k*a))
    end
    return phi
end
#at a distance x from the cylinder of radius r.
### BEM model

function TinyLev2D(r=8.575*5,k=40000*2*π/343000,l=10,h=40)
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    xc,yc = wavenurbs2D.calcula_centro(0,0,0,h,-r)
    POINTS =[1 0 0
       	     2 -l 0
       	     3 -l h 
       	     4 0 h
             5 2*xc 0
             6 2*xc h
             7 2*xc+l h
             8 2*xc+l 0
             ];
    SEGMENTS = [1 1 2 0
                2 2 3 0
                3 3 4 0
                4 4 1 r

                5 5 6 r              
                6 6 7 0
                7 7 8 0
                8 8 5 0
                ];
    BCSeg = [1 1 0 0
             2 1 0 0
             3 1 0 0
             4 1 1 0
             5 1 -1 0
             6 1 0 0
             7 1 0 0
             8 1 0 0];
    crv =  wavenurbs2D.format_dad_iso(POINTS,SEGMENTS)
    n_pint = 50; # Number of domain points
    PONTOS_dom = zeros(n_pint*n_pint,3);
    iniciox = POINTS[1,2];
    finalx = POINTS[5,2]
    passox = (finalx-iniciox)/n_pint
    inicioy = POINTS[1,3];
    finaly = POINTS[4,3];
    passoy = (finaly-inicioy)/n_pint
    global iter = 1
    for j = 1:n_pint
        for i = 1:n_pint
            PONTOS_dom[iter,:] = [iter iniciox+(i-1)*passox  inicioy+(j-1)*passoy];
            iter +=1
        end
    end
    fc = [0 0 0]
    info = [POINTS,SEGMENTS,BCSeg,k];
    t = @elapsed phi,q,phi_dom = wavenurbs2D.solveH(info, PONTOS_dom, fc, k)  
    #return norm(phi_dom.^2 .- phi_cylinder(k,r,PONTOS_dom[:,2]).^2)./size(PONTOS_dom,1)
    return phi,q,phi_dom, PONTOS_dom, crv
end
