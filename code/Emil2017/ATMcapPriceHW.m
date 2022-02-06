function price=ATMcapPriceHW(tenor,kappa,sigma,tenorZeroRates,zeroRates)
% # of prices 
N=size(tenor,1);

% Preallocate price(s)
price=nan(N,1);

% Calculate discount bonds
tenorZCB=(0.5:0.5:max(tenor))';
y=hermiteInterpolation(tenorZeroRates,zeroRates,tenorZCB);
ZCB=exp(-tenorZCB.*y);

% Preallocate par swap rates
parSwapRate=nan(N,1);
for i=1:N
    T=tenor(i);
    B=tenorZCB==T;
    numerator=1-ZCB(B);
    B=tenorZCB<=T&mod(tenorZCB,1)==0;
    denominator=sum(ZCB(B));
    parSwapRate(i)=numerator/denominator;
end

% Loop over tenors
for i=1:N
    % Calculate pseudo-strike
    K=1/(1+0.5*parSwapRate(i));
    
    % Calculate each caplet
    caplet=nan(tenor(i)*2,1);
    Bt=1;
    T=0;
    for j=1:size(caplet,1)
        Bs=ZCB(j);
        S=tenorZCB(j);
        caplet(j)=bondOptionHullWhite(Bs,Bt,K,0,T,S,kappa,sigma,'put');
        Bt=Bs;
        T=S;
    end
    price(i)=1/K*sum(caplet);
end

end