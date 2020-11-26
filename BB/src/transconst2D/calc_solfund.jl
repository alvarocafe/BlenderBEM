function  [G,H]=calc_solfund(RA,ISTEP,AT,CS)
#Calcula as soluções fundamentais

# THIS SUBROUTINE COMPUTES THE FUNDAMENTAL SOLUTION

G=0.0;
H=0.0;
A0=CS*ISTEP*AT/RA;
if(A0<=1.0)
    return
end
A1=CS*(ISTEP-1)*AT/RA;
A2=CS*(ISTEP-2)*AT/RA;
if(A1<=1.0)
    G=log(A0+raiz(A0))/(2*pi);
    H=-raiz(A0)/(2*pi)/CS/AT;
elseif(A2<=1.0)
    G=log((A0+raiz(A0))/(A1+raiz(A1)))/(2*pi);
    H=-((2-ISTEP)*(3*ISTEP-2)*(CS*AT/RA)^2+3.0)/(raiz(A0)+ 2.0*raiz(A1))/(2*pi)/CS/AT;
else
    G=log((A0+raiz(A0))/(A1+raiz(A1)))/(2*pi);
    H=-(2.0-2.0*(A0^2+A2^2-1.0)/(raiz(A0)*raiz(A2)+A0*A2))/(raiz(A0)+ 2.0*raiz(A1)+raiz(A2))/(2*pi)/CS/AT;
end
