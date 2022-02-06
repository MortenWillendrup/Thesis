%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Bootstrapping the borrower distribution
%
%   Emil Rode 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear
clc;clear;close all;

% Seed
seed=12345;
rng(seed)

% Time
T=120;

% Borrowers
N=10000;

% Some KNOWN function
t=(1:T)'./12;
Xstar=max(-0.0+cumsum(sqrt(0.001)*randn(size(t,1),1)),0);

% Number of bootstrap samples
B=1000;
THETA=nan(B,3);
wb=waitbar(0);
for h=1:B
    waitbar(h/B,wb)
    % Parameters A and B to be estimated
    Alpha=2;
    Beta=5;
    X=betarnd(Alpha,Beta,N,1);

    % Prepayment intensity to be estimated per time dt
    lambda=0.5;
    dt=0.25;

    % Prepayment indicator
    tau=false(N,T);
    U=rand(N,T);

    % For each borrower (note: independent!)
    for j=1:N
        % For each time unit
        for i=1:T
            % If this is fullfiled ...
            if X(j)<Xstar(i)
                % ... then prepayment happens with prob 1-exp(-lambda)
                if U(j,i)<1-exp(-lambda*dt)
                    tau(j,i)=true;
                    break;
                end
            end
        end
    end

    % Calculate prepayments
    prepayments=sum(tau,1)';
    Nt=N-cumsum(prepayments);
    prepayRates=prepayments./Nt;

    prepayProbs=nan(size(prepayRates));
    for i=1:size(Xstar,1)
        prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,Alpha,Beta,1/dt);
    end


    % Estimate
    theta=[Alpha;Beta;lambda];
    lb=zeros(3,1);
    ub=[10;10;3];
    THETA(h,:)=fmincon(@(theta)GMMform2(theta,prepayRates,Xstar),theta,[],[],[],[],lb,ub)';
end


XXX=tau-repmat(prepayProbs',N,1);
YYY=zeros(T,T);
for i=1:N
    YYY=YYY+XXX(i,:)'*XXX(i,:);
end
S=1/N*YYY;

% Calc D_T by finite difference
prepayProbsUp=nan(size(prepayRates));
prepayProbsDown=nan(size(prepayRates));
eps=0.00001;
for i=1:size(Xstar,1)
    prepayProbsUp(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,Alpha+eps,Beta,1/dt);
    prepayProbsDown(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,Alpha-eps,Beta,1/dt);
end
deltaA=(prepayProbsUp-prepayProbsDown)./(2*eps);

prepayProbsUp=nan(size(prepayRates));
prepayProbsDown=nan(size(prepayRates));
eps=0.00001;
for i=1:size(Xstar,1)
    prepayProbsUp(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,Alpha,Beta+eps,1/dt);
    prepayProbsDown(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,Alpha,Beta-eps,1/dt);
end
deltaB=(prepayProbsUp-prepayProbsDown)./(2*eps);

prepayProbsUp=nan(size(prepayRates));
prepayProbsDown=nan(size(prepayRates));
eps=0.00001;
for i=1:size(Xstar,1)
    prepayProbsUp(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda+eps,Alpha,Beta,1/dt);
    prepayProbsDown(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda-eps,Alpha,Beta,1/dt);
end
deltaLambda=(prepayProbsUp-prepayProbsDown)./(2*eps);


DT=[deltaA,deltaB,deltaLambda];

V=inv((DT'*inv(S)*DT));

Variance=1/N*V;

theta+sqrt(diag(Variance))*[norminv(0.025) norminv(0.975)]

subplot(1,3,1);histogram(THETA(1:h,1));
subplot(1,3,2);histogram(THETA(1:h,2));
subplot(1,3,3);histogram(THETA(1:h,3));

function [val,expectedPrepayRates]=GMMform2(theta,prepayments,Xstar)
% Get parameters
A=theta(1);
B=theta(2);
lambda=theta(3);
dt=0.25;

% Calculate prepayment probabilities
prepayProbs=nan(size(prepayments));
for i=1:size(Xstar,1)
    prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,A,B,1/dt);
end

% Calculate prepayment rates
remaining=1-cumsum(prepayProbs);
expectedPrepayRates=prepayProbs./remaining;

% Calculate quadratic form
f=prepayments-expectedPrepayRates;
val=f'*f;
end