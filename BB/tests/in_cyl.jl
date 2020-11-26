using Plots
pyplot()
include("../BEM_base.jl")
mshinfo = const3D_tri.lermsh("cylinder3D_tri.msh",3) #Read the mesh generated
BCFace = [1 1 0
          2 1 0
          3 0 1
          4 1 0];
# Build the domain points
L = 10; # Length of the vocal tract
n_pint = 40; # Number of domain points
#n_pint = 10; # Number of domain points
PONTOS_int = zeros(n_pint,4);
delta = 1; # distance from both ends 
passo = (L-2*delta)/(n_pint-1);
for i = 1:n_pint
    PONTOS_int[i,:] = [i 0 0 delta+(i-1)*passo];
end
k=1*pi/L

@elapsed uH,qH,uintH,qintH = const3D_tri.solveH(mshinfo,PONTOS_int,BCFace,k)
@elapsed u,q,uint,qint = const3D_tri.solve(mshinfo,PONTOS_int,BCFace,k)

plot(PONTOS_int[:,4],real(uint))
plot(PONTOS_int[:,4],real(uintH))
#pnts, elm, elmint, CDC = const3D_tri.lermsh("cylinder3D_tri.msh")
