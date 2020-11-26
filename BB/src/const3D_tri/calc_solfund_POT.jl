function calc_solfund_POT(x,y,z,xd,yd,zd,n,k)
rx=x-xd;
ry=y-yd;
rz=z-zd;

r =sqrt(rx^2+ry^2+rz^2);
Tast = 1.0/(4.0*pi*r);
qast = (rx*n[1] + ry*n[2] + rz*n[3])/(4.0*pi*r^3.0);

return Tast, qast
end

