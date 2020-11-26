function G=calcula_HeGs(CT0,CT1,SL)
# integração singular

#  THIS SUBROUTINE COMPUTES THE VALUE OF THE MATRIX G COEFFICIENTS
#  THAT RELATE AN ELEMENT WITH ITSELF


RS0=max(1.0,CT0/SL);
RS1=max(1.0,CT1/SL);
G=CT0*(acosh(RS0)/RS0+asin(1.0/RS0))/pi-CT1*(acosh(RS1)/RS1+ asin(1.0/RS1))/pi;
