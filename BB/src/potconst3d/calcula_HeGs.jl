function calcula_HeGs(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,k)
            
            
# Integra��o singular das matrizes H e G. No caso da matriz G, o elemento
# triangular � decomposto em tr�s quadril�teros degenerados na forma
# de tri�ngulos. Os dois primeiros n�s destes quadril�teros s�o
# coincidentes (formam um s� v�rtice do tri�ngulo) e correspondem ao ponto
# onde existe a singularidade, ou seja, ao ponto fonte que se localiza no
# centr�ide do elemento. Isto faz com que haja uma concentra��o de pontos
# de integra��o junto � singularidade, al�m do jacobiano ser igual a zero
# na singularidade. No caso da matriz H, a integra��o � anal�tica e sempre
# ser� igual a -1/2.

npg=length(qsi); # N�mero de pontos de integra��o
g = 0.0; # inicializa��o da matriz G

for kk=1:3
    x1t=xd; # coordenada x do primeiro n� do quadrilatero desgenerado
    y1t=yd; # coordenada y do primeiro n� do quadrilatero desgenerado
    z1t=zd; # coordenada z do primeiro n� do quadrilatero desgenerado
    x2t=xd; # coordenada x do segundo n� do quadrilatero desgenerado
    y2t=yd; # coordenada y do segundo n� do quadrilatero desgenerado
    z2t=zd; # coordenada z do segundo n� do quadrilatero desgenerado
    
    if(kk==1) # Terceiro e quarto n�s do primeiro quadril�tero desgenerado
        x3t=x1; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y1; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z1; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x2; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y2; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z2; # coordenada z do quarto n� do quadrilatero desgenerado
        
    elseif(kk==2) # Terceiro e quarto n�s do segundo quadril�tero 
                                  % desgenerado
        x3t=x2; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y2; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z2; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x3; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y3; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z3; # coordenada z do quarto n� do quadrilatero desgenerado
        
    elseif(kk==3) # Terceiro e quarto n�s do terceiro quadril�tero 
                                  % desgenerado
        x3t=x3; # coordenada x do terceiro n� do quadrilatero desgenerado
        y3t=y3; # coordenada y do terceiro n� do quadrilatero desgenerado
        z3t=z3; # coordenada z do terceiro n� do quadrilatero desgenerado
        x4t=x1; # coordenada x do quarto n� do quadrilatero desgenerado
        y4t=y1; # coordenada y do quarto n� do quadrilatero desgenerado
        z4t=z1; # coordenada z do quarto n� do quadrilatero desgenerado
        
    end
    
    for ii = 1: npg # la�o sobre a primeira vari�vel de integra��o
        for jj = 1: npg # la�o sobre a segunda vari�vel de integra��o
            [N1,N2,N3,N4] = calc_fforma(qsi[ii],qsi[jj]); % Fun��es de 
                                                       %  forma
            xc = x1t*N1+x2t*N2+x3t*N3+x4t*N4; # coordenada x do
            # ponto campo
            yc = y1t*N1+y2t*N2+y3t*N3+y4t*N4; # coordenada y do
            # ponto campo
            zc = z1t*N1+z2t*N2+z3t*N3+z4t*N4; # coordenada z do
            # ponto campo
            J = calc_jacobiano(x1t,y1t,z1t,x2t,y2t,z2t,x3t,y3t,z3t, ...
                x4t,y4t,z4t,qsi[ii],qsi[jj]); # jacobiano(varia ao longo 
                         #  do elemento desgenerado)
            Tast = calc_solfund(xd,yd,zd, xc, yc, zc, [0 0 0], k); # Sol.
            # fudamental de temperatura
            g = g + w[jj]*w[ii]*J*Tast; # integra��o num�rica da matriz G

        end
    end
end

h=-1/2 # Integra��o anal�tica da matriz H

return g,h

end













