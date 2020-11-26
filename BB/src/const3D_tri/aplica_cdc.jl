function aplica_cdc(G,H,CDC)
# Aplica as condi��es de contorno trocando as colunas das matrizes H e G
ncdc = length(CDC[:,1]); # n�mero de linhas da matriz CDC
A=copy(H);
B=copy(G);
for i=1:ncdc # La�o sobre as condi��es de contorno
    tipoCDC = CDC[i,2]; # Tipo da condi��o de contorno
    if tipoCDC == 0 # A temperatura � conhecida
        colunaA=-A[:,i]; # Coluna da matriz H que ser� trocada
        A[:,i]=-B[:,i]; # A matriz H recebe a coluna da matriz G
        B[:,i]=colunaA; # A mstriz G recebe a coluna da matriz H
    end
end;

#valoresconhecidos=complex(CDC[:,3],CDC[:,4]); # Valores das condicoes de contorno
valoresconhecidos = CDC[:,3];
b=B*valoresconhecidos; # vetor b

return A,b
end

