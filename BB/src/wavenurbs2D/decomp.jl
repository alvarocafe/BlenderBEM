function  bezierExtraction2D(uknot, vknot, p, q)
    Cxi, nb1, conn1, uu  = bezierExtraction(uknot, p);
    Cet, nb2, conn2, vv  = bezierExtraction(vknot, q);
  
    size1 = size(Cxi[:,:,1], 1);
    size2 = size(Cet[:,:,1], 1);
  
    C = zeros(size1 * size2, size1 * size2, nb1 * nb2);
    Conn = Array{Tuple}(nb1 * nb2)
    range = Array{Tuple}(nb1 * nb2, 2)
    for eta = 1:nb2
        for xi = 1:nb1
            e = (eta - 1) * nb1 + xi
            for row = 1:size2
                ird = (row - 1) * size1 + 1
                jrd =  row * size1
                for col = 1:size2
                    icd = (col - 1) * size1 + 1
                    jcd =  col * size1
                    C[ird:jrd,icd:jcd,e] = Cet[row,col,eta] * Cxi[:,:,xi];
          # display(C[ird:jrd,icd:jcd,e] )
                end
            end
            Conn[e] = (conn1[xi], conn2[eta])
            range[e,:] = [(uu[xi], uu[xi + 1]) (vv[eta], vv[eta + 1])]
      # for i =1:size1
      # for j =1:size2
      # conn[i+(j-1)*size1,e]=conn1[xi][i]+(conn2[eta][j]-1)*size1
      # end
      # end
        end
    end
    C, Conn, range
  
end


function bezierExtraction(knot, p)
    uk = unique(knot)
    nb1 = length(uk) - 1
    conn = Vector(undef, nb1)
    conn[1] = 1:p + 1
    for i = 2:nb1
        conn[i] = conn[i - 1] .+ sum(knot .== uk[i])
    end
    m  = length(knot) - p - 1;
    a  = p + 1;
    b  = a + 1;
    nb1 = length(uk) - 1
    C = zeros(p + 1, p + 1, nb1)
    C[:,:,1] = Matrix(1.0I, p + 1, p + 1)
    nb = 1
    while b <= m
        C[:,:,nb + 1] = Matrix(1.0I, p + 1, p + 1)
        i = b;
        while b <= m && knot[b + 1] == knot[b]
            b = b + 1;
        end
    
        multiplicity = b - i + 1;
        if multiplicity < p
            numerator = knot[b] - knot[a];
            alphas = zeros(p)
            for j = p:-1:multiplicity + 1
                alphas[j - multiplicity] = numerator / (knot[a + j] - knot[a]);
            end
            r = p - multiplicity;
            for j = 1:r
                save = r - j + 1;
                s = multiplicity + j;
                for k = p + 1:-1:s + 1
                    alpha = alphas[k - s];
                    C[:,k,nb] = alpha * C[:,k,nb] + (1 - alpha) * C[:,k - 1,nb];
                end
                if b <= m
                    C[save:save + j,save,nb + 1] = C[p - j + 1:p + 1,p + 1,nb]
                end
            end
            nb = nb + 1;
            if b <= m
                a = b;
                b = b + 1;
            end
        elseif multiplicity == p
            if b <= m
                                nb = nb + 1; a = b; b = b + 1;
            end
        end
    end
    C, nb, conn, uk
end


function bernsteinbasis(p, q, xi, eta)
  # Initialization
    ncpt = (p + 1) * (q + 1);
    B = zeros(ncpt, 1);
    dBdxi = zeros(ncpt, 2);
    for j = 1:q + 1
        for i = 1:p + 1
            B[(p + 1) * (j - 1) + i] = bernstein(p, i, xi) * bernstein(q, j, eta);
            dBdxi[(p + 1) * (j - 1) + i,1] = 0.5 * p * (bernstein(p - 1, i - 1, xi) - bernstein(p -
      1, i, xi)) * bernstein(q, j, eta);
            dBdxi[(p + 1) * (j - 1) + i,2] = bernstein(p, i, xi) * 0.5 * q * (bernstein(q - 1, j -
      1, eta) - bernstein(q - 1, j, eta));
        end
    end
    B, dBdxi
end

function bernstein(p, a, xi)
    if p == 0 && a == 1
        B = 1;
    elseif p == 0 && a != 1
        B = 0;
    else
        if a < 1 || a > p + 1
            B = 0;
        else
            B1 = bernstein(p - 1, a, xi);
            B2 = bernstein(p - 1, a - 1, xi);
            B = 0.5 * (1 - xi) * B1 + 0.5 * (1 + xi) * B2;
        end
    end
    B
end

function basisfundecomp2d(Be, dBedxi, Ce, we)
    Wb = Ce' * diag(we)
  # Bezier weight functions (denomenator of NURBS)
    wb        = dot(Be, Wb);
    dwbdxi = dot(dBedxi[:,1], Wb);
    dwbdxi = [dwbdxi dot(dBedxi[:,2], Wb)]
  # Shape function and derivatives
    R          = we * Ce * Be / wb;
    dRdxi = we * Ce * (dBedxi[:,1] / wb - dwbdxi[1] * Be / (wb * wb));
    dRdxi = [dRdxi we * Ce * (dBedxi[:,2] / wb - dwbdxi[2] * Be / (wb * wb))]
    R, dRdxi
end


function basisfundecomp(Be, dBedxi, Ce, w)
    Wb = Ce' * w
  # display(Be)
  # display(Wb)
  
    wb        = dot(Be, Wb);
    R         = diagm(0 => w) * Ce * Be / wb;
    dwbdxi = Vector(undef, 2)
    dwbdxi[1] = dot(dBedxi[:,1], Wb);
    dwbdxi[2] = dot(dBedxi[:,2], Wb);
  # Shape function and derivatives
    R          = diagm(0 => w) * Ce * Be / wb;
    dRdxi = diagm(0 => w) * Ce * (dBedxi[:,1] / wb - dwbdxi[1] * Be / (wb * wb));

    R, dRdxi
end
function evalsup(sup, u)
    uu = unique(sup.uknots)
    uv = unique(sup.vknots)
    indu = sum(u[1] .> uu)
    indv = sum(u[2] .> uv)
    ind = (length(uu) - 1) * (indv - 1) + indu
    unorm = 2 * (u[1] - uu[indu]) / (uu[indu + 1] - uu[indu]) - 1
    vnorm = 2 * (u[2] - uv[indv]) / (uv[indv + 1] - uv[indv]) - 1
    shapes, derivs = bernsteinbasis(sup.order[1] - 1, sup.order[1] - 1, unorm, vnorm);
  # control= reshape(sup.coefs[:,sup.conn[ind][1],sup.conn[ind][2]],4,length(sup.conn[ind][1])*length(sup.conn[ind][2]))
    control = reshape(sup.coefs[:,sup.conn[ind][1],sup.conn[ind][2]], 4, length(sup.conn[ind][1]) * length(sup.conn[ind][2]))
  
  # @show u
  # @show ind
    we = control[4,:]
    ptse = control[1:3,:] ./ we'
    R, dR = basisfundecomp2d(shapes, derivs, sup.C[:,:,ind], diagm(0 => we))
    x = ptse * R; # Calcula a coordenada x do ponto de integra��o
  #@show norm(x)
    dxdqsi = ptse * dR
    dgamadqsi = norm(cross(dxdqsi[:,1], dxdqsi[:,2]))        
    n = cross(dxdqsi[:,1], dxdqsi[:,2]) / dgamadqsi; # Componentes do vetor normal
  # Rcheio=zeros(sup.number[1]*sup.number[2])
  # Rcheio[(ones(Int64,length(sup.conn[ind][2]))*(sup.conn[ind][1]-1)'*(sup.number[2]).+sup.conn[ind][2])[:]]=R
  # x,Rcheio,n   
    x, R ./ we, n  
end

function evalcrv(crv, u)
    uu = unique(crv.knots)
    ind = sum(u .> uu)
    unorm = 2 * (u - uu[ind]) / (uu[ind + 1] - uu[ind]) - 1

    shapes, derivs = bernsteinbasis(crv.order - 1, 0, unorm, 0);
  # control= reshape(sup.coefs[:,sup.conn[ind][1],sup.conn[ind][2]],4,length(sup.conn[ind][1])*length(sup.conn[ind][2]))
    control = crv.coefs[:,crv.conn[ind]]
  
  # @show u
  # @show ind
    we = control[4,:]
    ptse = control[1:3,:] ./ we'
    R, dR = basisfundecomp2d(shapes, derivs, crv.C[:,:,ind], diagm(0 => we))
    x = ptse * R; # Calcula a coordenada x do ponto de integra��o
  #@show norm(x)
    dxdqsi = ptse * dR
    dgamadqsi = norm(dxdqsi);

  
    nx = dxdqsi[2] / dgamadqsi; # Componente x do vetor normal unit�rio
    ny = -dxdqsi[1] / dgamadqsi; # Componente y do vetor normal unit�rio
  # Rcheio=zeros(sup.number[1]*sup.number[2])
  # Rcheio[(ones(Int64,length(sup.conn[ind][2]))*(sup.conn[ind][1]-1)'*(sup.number[2]).+sup.conn[ind][2])[:]]=R
  # x,Rcheio,n   
  # x, R ./ we, [nx,ny]  
  x, R, [nx,ny]  
end