function [zero,forward,zeroTenors]=swapCalibration(swapTenors,swap,ciborTenors,cibor)
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

% Remove cibor >6M
B=ciborTenors>0.5;
ciborTenors(B)=[];
cibor(B)=[];

% Create zero tenors
zeroTenors=[ciborTenors;swapTenors];

% CIBOR is straight forward y=-1/tau*log(1/(1+delta*L))
ciborZero=-1./ciborTenors.*log(1./(1+ciborTenors.*cibor));

% Get number of swap tenors
N=size(swapTenors,1);

% Initial guess is just set to the swap rates
guess=swap;
swapFit=swap+1;

% Find all relevant tenors
tenors=(0.5:0.5:max(swapTenors))';
while sum(abs(swap-swapFit))*10000>1 % Max 0.1 BPS deviation
    % Interpolate zero rates
    zero=hermiteInterpolationFwd(zeroTenors,[ciborZero;guess],tenors);
    
    % Calculate ZCBs
    ZCB=exp(-tenors.*zero);

    % Calculate new zero rates through Hagan & West algorithm 
    sumP=nan(N,1);
    for i=1:N
        B=tenors<swapTenors(i);
        sumP(i)=sum(ZCB(B));
    end
    temp=(1-swap.*sumP*0.5)./(1+swap*0.5);
    
    % Compute new guess
    guess=-1./swapTenors.*log(temp);
    
    % Calculate theoretical swap rates
    swapFit=nan(size(swap,1),1);
    for i=1:N
        T=swapTenors(i);
        B=tenors==T;
        numerator=1-ZCB(B);
        B=ismember(tenors,0.5:0.5:T);
        denominator=0.5*sum(ZCB(B));
        swapFit(i)=numerator/denominator;
    end
end

% Compute output
[zero,forward]=hermiteInterpolationFwd(zeroTenors,[ciborZero;guess],zeroTenors);

end