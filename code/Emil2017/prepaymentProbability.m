function val=prepaymentProbability(x,Xstar,lambda,A,B,n,varargin)
% Sort data    
Xbar=sort([0;Xstar;1]);

% Preallocate
N=size(Xbar,1)-1;

% Prepayment probability
P=1-exp(-1/n*lambda);

% Calculate distribution 
summation=(1-P).^(N-(1:N)').*max(betacdf(min(x,Xbar(2:end)),A,B)-betacdf(Xbar(1:end-1),A,B),0);
%summation=(1-P).^(N-(1:N)').*max(betacdf(min(x,Xbar(2:end)),A,B)-betacdf(Xbar(1:end-1),A,B),0);

% Calculate output
val=P(end)*sum(summation);
end
