function calc_dfforma(qsi, eta)

dNdqsi = [1 0 -1]

dNdeta = [0 1 -1]
return dNdqsi, dNdeta
end


function calc_dfforma_quad(qsi, eta)

dNdqsi = (1/4)*[-(1 - eta);
                  (1 - eta);
                  (1 + eta);
                 -(1 + eta)];

dNdeta = (1/4)*[-(1 - qsi);
                 -(1 + qsi);
                  (1 + qsi);
                  (1 - qsi)];
return dNdqsi, dNdeta
end

