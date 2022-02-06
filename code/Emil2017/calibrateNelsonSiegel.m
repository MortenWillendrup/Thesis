function [b0,b1,b2,b3,t1,t2]=calibrateNelsonSiegel(ciborTenors,cibor,swapTenors,swap)
% Sort away nan
B=isnan(cibor);
cibor(B)=[];
ciborTenors(B)=[];
B=isnan(swap);
swap(B)=[];
swapTenors(B)=[];

% Preallocate parameters
theta=ones(6,1);

% Minimise squared errors
theta=fminunc(@(theta)squaredErrors(theta,ciborTenors,cibor,swapTenors,swap),theta);

% Extract results
b0=theta(1);
b1=theta(2);
b2=theta(3);
b3=theta(4);
t1=theta(5);
t2=theta(6);
end

function val=squaredErrors(theta,ciborTenors,cibor,swapTenors,swap)
% Create column data
ciborTenors=ciborTenors(:);
cibor=cibor(:);
swapTenors=swapTenors(:);
swap=swap(:);

% Calculate Zero Coupons for all tenors
tenors=(1:1:max(swapTenors))';
ZCB=nelsonSiegel(tenors,theta(1),theta(2),theta(3),theta(4),theta(5),theta(6));

% Fit CIBOR rates
P=nelsonSiegel(ciborTenors,theta(1),theta(2),theta(3),theta(4),theta(5),theta(6));
ciborFit=1./ciborTenors.*(1./P-1);

% Fit swap rates
P=nelsonSiegel(swapTenors,theta(1),theta(2),theta(3),theta(4),theta(5),theta(6));

swapFit=nan(size(swap,1),1);
for i=1:size(swap,1)
    T=swapTenors(i);
    numerator=1-P(i);
    B=tenors<=T;
    denominator=sum(ZCB(B));
    swapFit(i)=numerator/denominator;
end

% Calculate sum of squared errors
val=(sum((ciborFit-cibor).^2)*0+sum((swapFit-swap).^2)*100000);
end