function monta_Teq(CDC,x)
# Separa fluxo e temperatura

# ncdc = n�mero de linhas da matriz CDC
# T = vetor que cont�m as temperaturas nos n�s
# q = vetor que cont�m o fluxo nos n�s

ncdc = length(CDC[:,1]);
nnos = ncdc
T = complex(zeros(nnos,1))
q = complex(zeros(nnos,1))
for i=1:ncdc # La�o sobre as condi��es de contorno
    tipoCDC=CDC[i,2]; # Tipo da condi��o de contorno
    #valorCDC=complex(CDC[i,3],CDC[i,4]); # Valor da condi��o de contorno
		valorCDC = CDC[i,3];
		valorcalculado=x[i]; # Valor que antes era desconhecido
    if tipoCDC == 1 # Fluxo � conhecido
        T[i] = valorcalculado; # A temperatura � o valor calculado
        q[i] = valorCDC; # O fluxo � a condi�ao de contorno
    else # A temperatura � conhecida
        T[i] = valorCDC; # A temperatura � a condi�ao de contorno
        q[i] = valorcalculado; # O fluxo � o valor calculado
    end
end

return T,q
end

