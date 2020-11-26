function calc_dsolfund_POT(x,y,z,xd,yd,zd,n,k)
# Evaluates the fundamental solutions for the derivative of the Helmholtz equation.
# Determine the distance between source and field points
rx=x-xd;
ry=y-yd;
rz=z-zd;
r =sqrt(rx^2+ry^2+rz^2);

drdn=(rx*n[1] + ry*n[2] + rz*n[3])/r; # Unit normal derivative of the distance
phiast = (1/(4*pi*r^2))*drdn;
qastx = (1/(4*pi*r^3))*(n[1] - rx*drdn);
qasty = (1/(4*pi*r^3))*(n[2] - ry*drdn);
qastz = (1/(4*pi*r^3))*(n[3] - rz*drdn);
qast = sqrt(qastx^2 +qasty^2 +qastz^2)
return phiast,qast
end

