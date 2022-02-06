function [kappa,sigma,deviation]=capCalibration(tenor,quotes,tenorZeroRates,zeroRates)
% Force column input
tenor=tenor(:);
quotes=quotes(:);

B=isnan(quotes);
quotes(B)=[];
tenor(B)=[];

% Initial guess
theta=[0.001;0.005];

% Optimisation
lb=zeros(2,1);
ub=[2;10];
theta=fmincon(@(theta)minFun(theta,tenor,quotes,tenorZeroRates,zeroRates),theta,[],[],[],[],lb,ub);

% Unpack variables
sigma=theta(1);
kappa=theta(2);

% Ouput largest deviation
if nargout>2
    deviation=max(abs(quotes-ATMcapPriceHW(tenor,kappa,sigma,tenorZeroRates,zeroRates)*10000)./quotes);
end

end

function output=minFun(theta,tenor,quotes,tenorZeroRates,zeroRates)
% Unpack variables
sigma=theta(1);
kappa=theta(2);

% Get theoretical prices in basis points
prices=ATMcapPriceHW(tenor,kappa,sigma,tenorZeroRates,zeroRates)*10000;

% Compute squared errors
output=sum((prices-quotes).^2)*100000;
end