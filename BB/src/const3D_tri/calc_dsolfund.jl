function calc_dsolfund(x,y,z,xd,yd,zd,n,k)
# Evaluates the fundamental solutions for the derivative of the Helmholtz equation.

# Determine the distance between source and field points
rx=x-xd;
ry=y-yd;
rz=z-zd;
r =sqrt(rx^2+ry^2+rz^2);

kr = complex(0,-k*r);
drdn=(rx*n[1] + ry*n[2] + rz*n[3]); # Unit normal derivative of the distance
U=exp(kr)./complex(4*pi*r,0.);
phiast_x = (kr-complex(1.,0.)).*U.*complex(rx./r.^2,0.);
phiast_y = (kr-complex(1.,0.)).*U.*complex(ry./r.^2,0.);
phiast_z = (kr-complex(1.,0.)).*U.*complex(rz./r.^2,0.);
qast_x = (1. /(4*pi))*(-n[1]*(1. /r.^3+complex(0,1)*k./(r.^2)) + (3. /r.^3+3*complex(0,1)*k./(r.^2)-k^2. /(r)).*rx.*drdn./r.^2).*exp(-complex(0,1)*k*r);
qast_y = (1. /(4*pi))*(-n[2]*(1. /r.^3+complex(0,1)*k./(r.^2)) + (3. /r.^3+3*complex(0,1)*k./(r.^2)-k^2. /(r)).*ry.*drdn./r.^2).*exp(-complex(0,1)*k*r);
qast_z = (1. /(4*pi))*(-n[3]*(1. /r.^3+complex(0,1)*k./(r.^2)) + (3. /r.^3+3*complex(0,1)*k./(r.^2)-k^2. /(r)).*rz.*drdn./r.^2).*exp(-complex(0,1)*k*r);
return phiast_x,phiast_y,phiast_z,qast_x,qast_y,qast_z
end

