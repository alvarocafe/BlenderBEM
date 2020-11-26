## Test cases for potential problems
# This script tests the BEM models against known problems
#for which analytical solutions are readily obtained.
#################### Test case n ####################
### Problem description
### Analytical solution
### BEM model
### return error

#################### Test case 1 ####################
### Heat conduction on a square plate
# Let a square domain of length L [m] and heat
#conductivity k [W/(m.K)], initially at uniform constant
#temperature be subjected to a temperature gradient on
#its boundaries on opposites sides in the x direction.
### Analytical solution
# A linear solution solves the problem.
correction=1;
T_square(x,L) = 1 - x./L;
q_square(x,L) = -correction./L;
k_res(L) = n*pi/L;

### BEM model
function square(ne=10,L=1,k=1)
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
    #ne = 10;
    MESH = [1 ne
	    2 ne
	    3 ne
	    4 ne];
    BCSeg = [1 1 0
	     2 0 0
	     3 1 0
	     4 0 1];
    PONTOS_dom = [1 L/2 L/2]
    ## Formatting geometry for applying the method
    # Conventional method with full influence matricces
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg) # Apply the discretization technique and builds the problems matrices for the geometrical points, physical nodes, elements' connectivity and boundary conditions
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    tic()
    #Start----Problem---------------Method----------------post-processing
    #--------------------------------^You are here!----------------------
    T,q,T_dom=const2D.solvepot(info,PONTOS_dom,[0],BCSeg,k)
    t = toq()
    ϵ = abs.(sqrt.(((T_dom .- T_square(L/2,L)).^2)./T_square(L/2,L).^2))
    ## CBIE - Conventional Boundary Integral Equation
    collocCoord,nnos2,crv,dcrv,E,CDC = nurbs2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    info = [collocCoord,nnos2,crv,dcrv,E]
    tic()
    Tiso, qiso, T_domiso = nurbs2D.solvepot(info,PONTOS_dom,[0],CDC,k)
    tiso = toq()
    ϵiso = abs.(sqrt.(((T_domiso .- T_square(L/2,L)).^2)./T_square(L/2,L).^2))
    return t, ϵ, T, q, T_dom, tiso, ϵiso, Tiso, qiso, T_domiso
end
