# Boundary element method implementation for the Helmholtz equation using constant tridimensional triangular elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Modules necessary: SpecialFunctions.jl
include("../../src/const3D_tri/dep.jl") # Includes the dependencies
function const3D_tri(info)
	NOS_GEO,ELEM,elemint,CDC = info
	FR = 300*2*pi;
	CW = 343;
	k = FR/CW;
	# Gaussian quadrature - generation of points and weights [-1,1]
	npg=6; # Number of integration points
	qsi,w = Gauss_Legendre(0,1,npg) # Generation of the points and weights
	NOS = mostra_geoTRI(NOS_GEO,ELEM) #Generate the physical nodes for constant elements
	nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
	CCFace = ones(28,3);
	for i = 1:28
		CCFace[i,:] = [i 1 1]
	end
	CDC = gera_CDC(ELEM,CCFace); #Monta a matriz de condicoes de contorno
	n_pint = 100
	PONTOS_int = zeros(n_pint,4)
	for i = 1:n_pint
		PONTOS_int[i,:] = [i 0 0 (140/n_pint)*i]
	end
	println("Building G and H matrices...")
	@time G,H,phi_inc = cal_GeH(NOS,NOS_GEO,ELEM,k,qsi,w,0) #Compute the G and H matrices
	println("Applying boundary conditions to build A and b for the linear system...")
	@time A,b = aplica_cdc(G,H,CDC) #Applies the boundary conditions and returns matrix A and vector b for the linear system
	println("Solving the linear system...")
	@time x = A\b # Solves the linear system
	println("Separating acoustic pressure from flux...")
	@time phi,q = monta_Teq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux
	println("Solving for domain points.")
	T_pint=calc_T_pint(PONTOS_int,NOS_GEO,ELEM,phi,q,k,qsi,w,0)
end
