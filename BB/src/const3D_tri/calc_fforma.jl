function calc_fforma(qsi, eta)

#--------------------------------------------------------------------------
# Dados de entrada:
# qsi, eta - pontos de Gauss onde as fun��es de forma s�o calculadas.
# Dados de sa�da:
# [N] - Fun��es de forma para um elemento quadrilateral linear calculadas
# em (qsi, eta).
#--------------------------------------------------------------------------

N = (1/4.)*[(1. - qsi).*(1. - eta);
               (1. + qsi).*(1. - eta);
               (1. + qsi).*(1. + eta);
               (1. - qsi).*(1. + eta)];
N1=N[1];N2=N[2];N3=N[3];N4=N[4];
return N1,N2,N3,N4
end

