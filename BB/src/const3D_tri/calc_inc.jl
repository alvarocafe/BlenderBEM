function calc_inc(xd,yd,zd,FR,CW,inc)
#Funcao para calcular o valor da influencia de uma onda incidente no
#elemento estudado
k = complex(0,FR/CW); #Numero de onda
d = [inc[2,1] inc[3,1] inc[4,1]]; #direcao de propagacao da onda incidente /d/ = 1
p = [xd yd zd];
dp = dot(d,p);
A_inc = inc[5,1]; #amplitude da onda
phi_inc = A_inc*exp(k*dp);
return phi_inc
end

