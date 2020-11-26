# Boundary element method implementation for the Helmholtz and Laplace
#equations using constant  bidimensional elements
# Author: Álvaro Campos Ferreira - alvaro.campos.ferreira@gmail.com
# Contains the dependencies for the linear constant element integration.
#The main function is const2D.solve() which builds the influence matrices,
#applies the boundary conditions, solves the linear system and returns the
#value of the potential and its gradient at boundary and domain points.
module wavenurbs2D

# include("Potencial.jl")
include("dad.jl")
include("decomp.jl")
include("nurbs.jl")
include("CalcHeG.jl")
include("telles.jl")
include("arvore.jl")
include("formatiso.jl")
include("beminterp.jl")
using FastGaussQuadrature, LinearAlgebra, SparseArrays, Statistics, KrylovMethods, Test, SpecialFunctions, Plots
pyplot()

function solveH(info, PONTOS_dom, fc, kmat)
    PONTOS,SEGMENTOS,CCSeg,k = info
    crv=format_dad_iso(PONTOS,SEGMENTOS)# formata os dados
    n = length(crv);	# N�mero total de elementos
    p=0;#refinamento p
    for i=1:n
        degree=crv[i].order-1
        coefs,knots = bspdegelev(degree,crv[i].coefs,crv[i].knots,p)
        crv[i] = nrbmak(coefs,knots)
    end
    h=5;#refinamento h
    for i=1:n
        novosnos=range(0,stop=1,length=h+2)
        degree=crv[i].order-1
        coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
        crv[i] = nrbmak(coefs,knots)
    end
    Tree1,Tree2,block= cluster(crv,max_elem=8,η = 1.0)#cluster(crv, max_elem=3,η = 1.0)
    indfonte,indcoluna,indbezier,tipoCDC,valorCDC,E,collocCoord,collocPts=indices(crv,CCSeg)
    HA,bi=Hinterp(indfonte,indbezier,indcoluna,E,Tree1,Tree2,block,crv,kmat,tipoCDC,valorCDC,collocCoord)
    xi,f = gmres(vet->matvec(HA,vet,block,Tree1,Tree2,indcoluna),bi,5,tol=1e-5,maxIter=1000,out=0) #GMRES nas matrizes hierarquicas
    Tc,qc=monta_Teq(tipoCDC,valorCDC,crv, xi) # Separa temperatura e fluxo
    T=E*Tc
    q=E*qc
    Tdom = calc_pintpot(PONTOS_dom,indcoluna,indbezier, crv,collocCoord, kmat,Tc,qc)
    return T,q,Tdom
end # end function solveH

function solve(info, PONTOS_dom, fc, kmat)
    PONTOS,SEGMENTOS,CCSeg,k = info
    crv=format_dad_iso(PONTOS,SEGMENTOS)# formata os dados
    n = length(crv);	# N�mero total de elementos
    p=0;#refinamento p
    for i=1:n
        degree=crv[i].order-1
        coefs,knots = bspdegelev(degree,crv[i].coefs,crv[i].knots,p)
        crv[i] = nrbmak(coefs,knots)
    end
    h=5;#refinamento h
    for i=1:n
        novosnos=range(0,stop=1,length=h+2)
        degree=crv[i].order-1
        coefs,knots = bspkntins(degree,crv[i].coefs,crv[i].knots,novosnos[2:end-1])
        crv[i] = nrbmak(coefs,knots)
    end
    indfonte,indcoluna,indbezier,tipoCDC,valorCDC,E,collocCoord,collocPts=indices(crv)

    A,B=CalcAeb(indfonte,indcoluna,indbezier, crv, kmat,E,tipoCDC)
    b=(B*(E\valorCDC))[:]
    x=A\b

    Tc,qc=monta_Teq(tipoCDC,valorCDC,crv, x) # Separa temperatura e fluxo
    T=E*Tc
    q=E*qc

    Tdom = calc_pintpot(PONTOS_dom,indcoluna,indbezier, crv, kmat,Tc,qc)

    return T,q,Tdom
end # end function solve
end # end module wavenurbs2D
