clc; clear;close all;
A=1;
B=4;
N=10000;
X=betarnd(A,B,N,1);
lambda=2;
n=4;
P=1-exp(-1/n*lambda);
U=rand(N,1)<P;
U2=rand(N,1)<P;
U3=rand(N,1)<P;
Xstar=[0.5;0.2;0.15];

tau1=X<Xstar(1)&U;
tau2=X<Xstar(2)&U2;
tau3=X<Xstar(3)&U3;
b=(~tau1)&(~tau2)&(~tau3);

temp=X(b);

% Prob prepay time 3
P1=sum(tau1)/N;
P2=sum(tau2&~tau1)/N;
P3=sum(tau3&~tau2&~tau1)/N;
Nt=[1,1-cumsum([P1 P2 P3])]';
1-Nt(2:end)./Nt(1:end-1)

prepayProbs=nan(3,1);
for i=1:size(Xstar,1)
    prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,A,B,n);
end
Nt=[1;1-cumsum(prepayProbs)];
prepayRates2=1-Nt(2:end)./Nt(1:end-1)


% Prob not prepay
sum(~tau3&~tau2&~tau1)/N;


prepaymentProbability(Xstar(3),Xstar(1:2),lambda,A,B,n)

x=123;
% 
% 
% histogram(X,200)
% hold on
% histogram(X(~tau1),200)
% histogram(X(~tau1&~tau2),200)
% histogram(X(~tau1&~tau2&~tau3),200)
% yyaxis right
% plot(xx,betapdf(xx,A,B))
% plot(xx(xx>0.3),betapdf(xx(xx>0.3),A,B),'linewidth',2)
% plot(xx(xx<0.3),(1-lambda)*betapdf(xx(xx<0.3),A,B),'linewidth',2)
% plot(xx(xx<0.2),(1-lambda)^2*betapdf(xx(xx<0.2),A,B),'linewidth',2)
% plot(xx(xx<0.15),(1-lambda)^3*betapdf(xx(xx<0.15),A,B),'linewidth',2)

% yyaxis right
% histogram(temp)
% yyaxis left
% 
PDF=@(x)(1-P)^3*betapdf(x,A,B).*(x<Xstar(3))...
    +(1-P)^2*betapdf(x,A,B).*(x<=Xstar(2)&x>=Xstar(3))...
    +(1-P)*betapdf(x,A,B).*(x<=Xstar(1)&x>Xstar(2))...
    +betapdf(x,A,B).*(x>Xstar(1));
% 
% plot(xx,PDF(xx))
% 
% 
% CDF=@(x)(1-lambda)^3*betacdf(x,A,B).*(x<Xstar(3))...
%     +(1-lambda)^2*betacdf(x,A,B).*(x<=Xstar(2)&x>=Xstar(3))...
%     +(1-lambda)*betacdf(x,A,B).*(x<=Xstar(1)&x>Xstar(2))...
%     +betacdf(x,A,B).*(x>Xstar(1));


Xbar=[0;sort(Xstar);1];
B2=P.^(size(Xstar,1):-1:0)';
XbarProb=betacdf(Xbar(1:end-1),A,B);
int=1/integral(PDF,0,1);
CDF=@(x)int*sum((betacdf(x,A,B)-XbarProb).*(x>=Xbar(1:end-1)&x<=Xbar(2:end)).*B2);

[f,x] = ecdf(temp);
plot(x,f,'linewidth',3)
hold on
%yy=(0:0.001:0.15)';
%plot(yy,int*(1-lambda)^3*betacdf(yy,A,B),'linewidth',3)
xx=(0.01:0.01:1)';
probs=nan(size(xx));
for i=1:size(probs,1)
    probs(i)=1/P*int*CDF2(xx(i),Xstar,lambda,A,B,n);
end
plot(xx,probs,'linewidth',2,'linestyle','--')


h=histogram(X(~tau1&~tau2&~tau3));
xxx=h.BinCounts;
plot(cumsum(h.BinCounts)./sum(~tau1&~tau2&~tau3));

histogram(CDF,100)






sum(~tau1)/N
P1=1-betacdf(Xstar(1),A,B)+(1-lambda)*betacdf(Xstar(1),A,B);
P2=1/P1*(1-betacdf(Xstar(2),A,B)+(1-lambda)*betacdf(Xstar(2),A,B));


PPP=betacdf(Xstar(3),A,B)-betacdf(min(Xstar(1:2)),A,B)+...
    betacdf(min(Xstar([1;3])),A,B)-betacdf(Xstar(2),A,B)*(1-lambda)+...
    betacdf(min(Xstar),A,B)*(1-lambda)^2
    

P=[(1-betacdf(Xstar,A,B)),(1-lambda)*betacdf(Xstar,A,B)];

P(1,1)*P(2,1)+P(1,1)*P(2,2)+P(1,2)*P(2,2)+P(1,2)*P(2,1)


SUM=P(1,1)*P(2,1)*P(3,1);
SUM=SUM+P(1,1)*P(2,1)*P(3,2);
SUM=SUM+P(1,1)*P(2,2)*P(3,2);
SUM=SUM+P(1,2)*P(2,2)*P(3,2);
SUM=SUM+P(1,2)*P(2,2)*P(3,1);
SUM=SUM+P(1,2)*P(2,1)*P(3,1)


