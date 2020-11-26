function calc_solfund(x,y,z,xd,yd,zd,n,k)
# Evaluates the fundamental solutions of the Helmholtz equation.

# Determine the distance between source and field points
rx=x-xd;
ry=y-yd;
rz=z-zd;
r =sqrt(rx^2+ry^2+rz^2);

ZW=complex(0.,k*r);
U=exp(ZW)/complex(4*pi*r,0.); # Fundamental solution for the velocity potential
drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r;	# Unit normal derivative of the distance
Q=(ZW-complex(1.,0.))*U*complex(drdn/r,0.); # Fundamental solution for the flux
return U,Q
end

