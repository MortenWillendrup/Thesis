%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Generate parameters to be used for estimation
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;

% Load market
mkt=market('date','2017-06-01');

% Extract rates
startDate='2002-01-01';
endDate='2017-06-01';
[dateVector,rates]=timeSeries(mkt,'swap1Y',startDate,endDate);

% Remove nan observations
B=isnan(rates);
dateVector(B)=[];
rates(B)=[];

% MLE
data=[dateVector,rates];
parameters=[0.3;0.05;0.01]; % Kappa, theta & sigma 
parameters=fmincon(@(parameters)MLE(parameters,data),parameters,[],[]);
kappa=parameters(1);

% Get implied vol
dateVector=(datenum(startDate):datenum(endDate))';
impliedVol=nan(size(dateVector));
for i=1:size(dateVector,1)
    mkt.date=dateVector(i);
    impliedVol(i)=capCalibrationKnownKappa(mkt,kappa);
end


% Save data
clear data
data.kappa=repmat(kappa,size(dateVector,1),1);
data.sigma=impliedVol;
data.dateVector=dateVector;
save('/users/helmig/dropbox/ku/speciale/Data/parameters','data')

function [logLik,errors]=MLE(parameters,data)
% Get parameters
kappa=parameters(1);
theta=parameters(2);
sig=parameters(3);

% Get time and rates
time=data(:,1);
rates=data(:,2);

% Generate likelihood
logLik=0;
errors=nan(size(data,1)-1,1);
for i=2:size(data,1)
    dt=(time(i)-time(i-1))/365;
    mu=exp(-kappa*dt)*rates(i-1)+theta*(1-exp(-kappa*dt));
    %nu=sig*exp(-kappa*dt)*sqrt(dt);
    nu=sig*sqrt(1-exp(-2*kappa*dt))/sqrt(2*kappa);
    logLik=logLik+...
        log(nu)+(rates(i)-mu)^2/(2*nu^2);
    
    errors(i-1)=(rates(i)-mu)/nu;
end
end