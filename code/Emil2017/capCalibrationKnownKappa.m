% Function that calibrates cap market
function sigma=capCalibrationKnownKappa(mkt,kappa)
% Initial guess
sigma=0.001;

% Optimisation
lb=0;
ub=10;
options=optimoptions('fmincon','display','none');
sigma=fmincon(@(sigma)minFun(sigma,kappa,mkt.capTenors,...
                             mkt.capCurve,(0.5:0.5:30)',...
                             mkt.zeroCurve((0.5:0.5:30)')),...
                             sigma,[],[],[],[],lb,ub,[],options);
end

% Sum of squared errors for cap calibration
function output=minFun(sigma,kappa,tenor,quotes,tenorZeroRates,zeroRates)
% Get theoretical prices in basis points
prices=ATMcapPriceHW(tenor,kappa,sigma,tenorZeroRates,zeroRates)*10000;

% Compute squared errors
output=sum((prices-quotes).^2)*100000;
end