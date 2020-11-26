function calc_jacobiano(x1,y1,z1,x2,y2,z2,x3,y3,z3,qsi,eta)
dNdqsi, dNdeta = calc_dfforma(qsi,eta); # Calcula a derivada das fun��es
  # de forma
dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3]
dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3]
dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3]

dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3]
dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3]
dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3]

g1 = dydqsi*dzdeta - dzdqsi*dydeta;
g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
g3 = dxdqsi*dydeta - dydqsi*dxdeta;
J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);

return J
end

function calc_jacobiano_quad(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,qsi,eta)
dNdqsi, dNdeta = calc_dfforma_quad(qsi,eta); # Calcula a derivada das fun��es
  # de forma
dxdqsi = x1*dNdqsi[1]+x2*dNdqsi[2]+x3*dNdqsi[3]+x4*dNdqsi[4];
dydqsi = y1*dNdqsi[1]+y2*dNdqsi[2]+y3*dNdqsi[3]+y4*dNdqsi[4];
dzdqsi = z1*dNdqsi[1]+z2*dNdqsi[2]+z3*dNdqsi[3]+z4*dNdqsi[4];

dxdeta = x1*dNdeta[1]+x2*dNdeta[2]+x3*dNdeta[3]+x4*dNdeta[4];
dydeta = y1*dNdeta[1]+y2*dNdeta[2]+y3*dNdeta[3]+y4*dNdeta[4];
dzdeta = z1*dNdeta[1]+z2*dNdeta[2]+z3*dNdeta[3]+z4*dNdeta[4];

g1 = dydqsi*dzdeta - dzdqsi*dydeta;
g2 = dzdqsi*dxdeta - dxdqsi*dzdeta;
g3 = dxdqsi*dydeta - dydqsi*dxdeta;
J = sqrt(g1^2.0 + g2^2.0 + g3^2.0);

return J
end

