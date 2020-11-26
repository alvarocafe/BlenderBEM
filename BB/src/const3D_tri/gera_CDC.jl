function gera_CDC(ELEM,CCFace)

# Gera a matriz de CDC
# CDC = [número do elemento, tipo da CDC, valor da CDC no nó 1,...
#                               valor da CDC no nó 2, valor da CDC no nó 3]
nelem = length(ELEM[:,1]);
CDC = zeros(nelem,3);
for i=ELEM[:,1]
    CDC[i,1] = i;
    CDC[i,2] = CCFace[ELEM[i,5],2];
    CDC[i,3] = CCFace[ELEM[i,5],3];
end
return CDC
end

