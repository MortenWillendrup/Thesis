function [zero,forward,zeroTenors]=swapCalibration2(swapTenors,swap,ciborTenors,cibor)
% Remove nan
B=isnan(cibor);
cibor(B)=[];
ciborTenors(B)=[];
B=isnan(swap);
swap(B)=[];
swapTenors(B)=[];

% Sort data
[ciborTenors,IX]=sort(ciborTenors);
cibor=cibor(IX);
[swapTenors,IX]=sort(swapTenors);
swap=swap(IX);

% Preallocate all tenors needed
tenors=[ciborTenors;(min(swapTenors):0.5:max(swapTenors))'];

% Initial guess
zero=[cibor;swap];
zeroTenors=[ciborTenors;swapTenors];

% Optimize quasi-newton
options=optimoptions('fminunc','algorithm','quasi-newton');
zero=fminunc(@(zero)minFun(zero,zeroTenors,swapTenors,swap,ciborTenors,cibor,tenors),zero,options);

if nargout>1
    [~,forward]=hermiteInterpolation(swapTenors,zero,swapTenors);
end
    
end

function val=minFun(zero,zeroTenors,swapTenors,swap,ciborTenors,cibor,tenors)
% All tenors from data

% Calculate interpolated zero rates for all needed tenors
zeroRates=hermiteInterpolation(zeroTenors,zero,tenors);

% Calculate zero coupon prices
ZCB=exp(-tenors.*zeroRates);

% Calculate theoretical swap rates
swapFit=nan(size(swap,1),1);
for i=1:size(swap,1)
    T=swapTenors(i);
    B=tenors==T;
    numerator=1-ZCB(B);
    B=(tenors<=T)&(tenors>=0.5);
    denominator=0.5*sum(ZCB(B));
    swapFit(i)=numerator/denominator;
end

% Calculate theoretical cibor rates
ciborFit=nan(size(cibor,1),1);
for i=1:size(cibor,1)
    ciborFit(i)=1/ciborTenors(i)*(1/ZCB(i)-1);
end

% Produce squared errors
val=(sum((cibor-ciborFit).^2)+sum((swap-swapFit).^2))*100000;
end