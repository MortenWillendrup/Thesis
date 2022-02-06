function zeroPrices=nelsonSiegel(mat,b0,b1,b2,b3,t1,t2)
exponent=b0*mat...
         +b1*t1*(1-exp(-1/t1*mat))...
         -b2*mat.*exp(-1/t1*mat)...
         +b2*t1*(1-exp(-1/t1.*mat))...
         -b3*mat.*exp(-1/t2*mat)...
         +b3*t2*(1-exp(-1/t2*mat));
zeroPrices=exp(-exponent);
end