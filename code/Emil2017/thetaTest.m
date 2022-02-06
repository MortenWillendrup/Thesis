

clc;clear;close all;
dt=0.5;
tenors=(dt:dt:30)';
mkt=market('date','2017-06-01');
curve=mkt.zeroCurve(tenors);
N=400;
r_min=-0.1;
r_max=0.5;
dr=(r_max-r_min)/N;
space=(r_min:dr:r_max)';

% Preallocate theta
theta=nan(size(tenors));

% Solve iteratively for theta
options=optimoptions('fminunc','Algorithm','quasi-newton');
kappa=mkt.kappa;
sigma=mkt.sigma;
shortRate=mkt.zeroCurve(1/12);
% for i=1:size(tenors,1)
%     if i==1
%         guess=0;
%     else
%         guess=theta(i-1);
%     end
%     marketPrice=exp(-tenors(i)*curve(i));
%     theta(i)=...
%         fminunc(@(guess)minFunc(guess,kappa,theta(1:i-1),sigma,space,dr,dt,N,i,marketPrice,shortRate),guess,options);
%     
% end


% Get forward curve and its derivative

forward=mkt.forwardCurve(mkt.marketTenors);
derivative=mkt.forwardCurveDiff(mkt.marketTenors);

% Calculate theta
theta2=forward+1/kappa*derivative...
    +sigma^2/(2*kappa)*(1-exp(-2*kappa*mkt.marketTenors));
theta2=hermiteInterpolation(mkt.marketTenors,theta2,tenors);

modelPrice=nan(size(tenors));
for i=1:size(tenors,1)
    modelPrices=priceZCB(kappa,theta2,sigma,space,dr,dt,N,i);
    modelPrice(i)=hermiteInterpolation(space,modelPrices,shortRate);
end

marketPrice=exp(-tenors.*curve);
plot([-1./tenors.*log(modelPrice),curve])
plot([modelPrice,exp(-tenors.*curve)])

function val=minFunc(guess,kappa,theta,sigma,space,dr,dt,N,J,marketPrice,shortRate)
    thetaNew=[theta;guess];
    modelPrices=priceZCB(kappa,thetaNew,sigma,space,dr,dt,N,J);
    modelPrice=hermiteInterpolation(space,modelPrices,shortRate);
    val=(modelPrice-marketPrice)^2*1000000;
end




function price=priceZCB(kappa,theta,sigma,space,dr,dt,N,J)
price=ones(N+1,1);
for j=J:-1:1
    % Mu and sigma
    Mu=kappa*(theta(j)-space);
    SIG=sigma^2;

    % A B C D vectors
    A=1/(4*dr)*Mu-1/(4*dr^2)*SIG;A=A(:);
    B=1/dt+1/2*SIG*1/(dr^2)+1/2*space;B=B(:);
    C=-1/(4*dr)*Mu-1/(4*dr^2)*SIG;C=C(:);
    D=1/dt-1/(2*dr^2)*SIG-1/2*space;D=D(:);

    % Correcting boundaries (zero convexity)
    B(1)=1/dt+Mu(1)/(2*dr)+1/2*space(1);
    C(1)=-Mu(1)/(2*dr);
    D(1)=1/dt-Mu(1)/(2*dr)-1/2*space(1);
    A(end)=Mu(end)/(2*dr);
    B(end)=1/dt-Mu(end)/(2*dr)+1/2*space(end);
    D(end)=1/dt+Mu(end)/(2*dr)-1/2*space(end);

    % Left Hand Side Matrix
    LHS1=[[zeros(1,N);diag(A(2:end))],zeros(N+1,1)];
    LHS2=diag(B);
    LHS3=[zeros(N+1,1),[diag(C(1:end-1));zeros(1,N)]];
    LHS=LHS1+LHS2+LHS3;

    % Right Hand Side Matrix
    RHS1=[[zeros(1,N);diag(-A(2:end))],zeros(N+1,1)];
    RHS2=diag(D);
    RHS3=[zeros(N+1,1),[diag(-C(1:end-1));zeros(1,N)]];
    RHS=RHS1+RHS2+RHS3;

    % Invert LHS
    invLHS=(inv(LHS));
    
    % Calculate new values of liabilities
    f=RHS*price;
    price=invLHS*f;
end
end