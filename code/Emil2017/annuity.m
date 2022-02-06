function schedule=annuity(coupon,principal,years,paymentsPerYear)
n=paymentsPerYear;
T=floor(years*n)/n;
%residual=years-T;
R_tilde=coupon/n;
Y=R_tilde/(1-(1+R_tilde)^(-n*T))*principal;
i=(0:(T*n))';
F_t=(1+R_tilde).^i*principal-Y*((1+R_tilde).^i-1)./R_tilde;
interest=[0;F_t(1:end-1)*R_tilde];
amortisation=[0;F_t(1:end-1)-F_t(2:end)];
payDates=(0:1/n:T)';%+residual;
payments=[0;repmat(Y,n*T,1)];
schedule=[payDates,interest,amortisation,payments,F_t];
end