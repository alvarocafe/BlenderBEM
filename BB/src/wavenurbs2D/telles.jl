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

eta = A*gamm.^3 .+ B*gamm.^2 .+ C*gamm .+ D;
Jt = 3*A*gamm.^2 .+ 2*B*gamm .+ C;
return  eta,Jt
end
