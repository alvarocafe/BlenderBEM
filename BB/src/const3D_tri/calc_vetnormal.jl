function calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3)
# Function que calcula o vetor unit�rio normal ao elemento

v1 = [x3,y3,z3] - [x2,y2,z2]; # vetor formado pela aresta 32 do elemento
v2 = [x1,y1,z1] - [x2,y2,z2]; # vetor formado pela aresta 12 do elemento
n = cross(v1, v2); # Produto vetorial entre v1 e v2 (vetor normal ao
                           # elemento)
#println(v1,v2,n)
if sum(abs.(n)) > 0.000001
	n = n./norm(n); # vetor unit�rio normal ao elemento
end
#println(v1,v2,n)
return n
end

