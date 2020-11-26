include("tests/wave_tests.jl")
using Plots
pyplot()
function contour_plt(n,k,r)
    pontos = zeros(n*n,3)
    jump = -5*r + π/16
    passox = r/10
    iter = 0
    for i in 1:n
        for j in 1:n
            iter += 1
            pontos[iter,:] = [iter jump+i*passox jump+j*passox]
        end
    end
    phi_dom = complex(zeros(n*n))
    [phi_dom[i] = phi_cylinder(k,r,sqrt(pontos[i,2]^2 + pontos[i,3]^2)) for i in 1:size(pontos,1)]
    return pontos, phi_dom
end

n = 100
k = 5
r = 0.5
pontos, phi_dom = contour_plt(n,k,r)

plot(pontos[:,2],pontos[:,3],real(phi_dom),st=:surface,c=:viridis)
#plot!(-pontos[:,2],pontos[:,3],real(phi_dom),st=:surface)
#plot!(-pontos[:,2],-pontos[:,3],real(phi_dom),st=:surface)
#plot!(pontos[:,2],-pontos[:,3],real(phi_dom),st=:surface,camera=(90,90))
