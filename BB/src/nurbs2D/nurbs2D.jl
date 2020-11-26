# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
#Start----Problem--------------Method-----------------post-processing
#-----------^You are here!-------------------------------------------

module nurbs2D
using SpecialFunctions, KrylovMethods
include("format.jl") # curve interpolation formatting
include("cal.jl") # element integration calculating functions
include("H_mat.jl") # H-Matrices support for building the cluster tree and blocks
include("interp.jl") # approximation  using Lagrange polynomial interpolation
include("ACA.jl") # approximation using ACA

function solve(info,PONTOS_int,fc,CDC,k)
    ## CBIE - Conventional Boundary Integral Equation
    collocCoord,nnos,crv,dcrv,E = info
    npg=6 # Number of integration points
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = cal_Aeb(1,1,collocCoord,nnos,crv,dcrv,E,k,CDC)
    #    A,b = cal_Aeb(b1,b1, [NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    x = A\b # Solves the linear system
    phi,qphi = monta_Teq(CDC,x)
    phi_dom = calc_phi_pint(PONTOS_int,collocCoord,nnos,crv,dcrv,k,phi,qphi)
    return phi, qphi, phi_dom, phi_dom
end

function solveH(info,PONTOS_int,fc,BCFace,k)
    ## H-Matrix BEM - Interpolation using Lagrange polynomial
    NOS,NOS_GEO,ELEM,CDC = info
    Tree,block = cluster(NOS[:,2:3],floor(sqrt(length(NOS))),2)
    # Gaussian quadrature - generation of points and weights [-1,1]
    npg=6 # Number of integration points
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    A,b = Hinterp(Tree,block,[NOS,NOS_GEO,ELEM,fc,qsi,w,CDC,k])
    x = gmres(vet->matvec(Ai,vet,block,Tree),bi,5,tol=1e-5,maxIter=1000,out=0)
    phi,qphi = monta_phieq(CDC,xi[1]) # Applies the boundary conditions
    phi_pinti = calc_phi_pint(PONTOS_int,NOS_GEO,ELEM,phii,qphii,fc,fc,qsi,w,k)
    return phii,qphii,phi_pinti,phi_pinti
end

function solvepot(info,PONTOS_int,fc,CDC,k)
    ## CBIE - Conventional Boundary Integral Equation
    collocCoord,nnos,crv,dcrv,E = info
    npg=6 # Number of integration points
    qsi,w = Gauss_Legendre(-1,1,npg) # Generation of the points and weights
    # H,G = nurbs2D.calcula_iso_POT(collocCoord,nnos,crv,dcrv,E,k)
    # A,b = aplica_CDC(G,H,CDC,E)
    A,b = cal_Aebpot(1,1,collocCoord,nnos,crv,dcrv,E,k,CDC)
    x = A\b # Solves the linear system
    T,q = monta_Teq(CDC,x)
    T_dom = calc_phi_pintpot(PONTOS_int,collocCoord,nnos,crv,dcrv,k,T,q)
    return T, q, T_dom, T_dom
end

end
