%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Pool example
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;

% Set parameters
kappa=0.29368;
mu=0.07935;
sig=0.11425;
q=-0.12165;
F0=100;
R=0.125;
n=4;
T=30;
lambda1=0;
lambda2=0.6452;

% Payment schedule
schedule=annuity(R,F0,T,n);
time=schedule(:,1);
F=schedule(:,end);

% Simulate short rate
r=nan(size(time));
r(1)=mu;
dt=time(2);
W=cumsum(randn(size(time,1),1)*sqrt(dt));
for i=2:size(time,1)
    r(i)=r(i-1)+kappa*(mu-r(i-1))*dt+sig*sqrt(r(i-1))*(W(i)-W(i-1));
end
% r(1:end)=0.1;
% r(21:end)=0.05;
% r(31:end)=0.1;
% r(41:end)=0.04;
% r(51:end)=0.1;
% r(61:end)=0.03;
% r(71:end)=0.1;
% r(81:end)=0.02;
% r(91:end)=0.1;
% r(101:end)=0.01;

% Pool
N=1000;
A=0.9618;
B=4.2268;
X=betainv(rand(N,1),A,B);
AO=nan(size(time));
AO(1)=N*F0;
tau=zeros(N,size(time,1));

% Initiate model
model=cStanton('kappa',kappa,...
               'mu',mu,...
               'sigma',sig,...
               'q',q,...
               'F0',F0,...
               'R',R,...
               'n',n,...
               'T',T,...
               'lambda1',lambda1,...
               'lambda2',lambda2);

% Generate prepayments
U=rand(size(tau));
h=waitbar(0);
optimality=zeros(N,size(time,1));
for j=1:N
    waitbar(j/N,h);
    model.X=X(j);
    [~,~,space,t,asset,liab]=model.pricingFD;
    if j==1
        [~,index1]=min(abs(space-r'));
        index2=nan(size(time));
        for i=1:size(time,1)
            [~,index2(i)]=min(abs(t-time(i)));
        end
    end
    
    for i=1:size(time,1)
        if liab(index1(i),index2(i))>=F(i)*(1+X(j))
            tau(j,i)=U(j,i)<1-exp(-dt*lambda2);
            if tau(j,i)==1
                break
            end
        end
    end
    for i=1:size(time,1)
        if liab(index1(i),index2(i))>=F(i)*(1+X(j))
            optimality(j,i)=1;
        end
    end
end
close(h)
tau2=cumsum(tau,2);
AO=nan(size(time));
A=nan(size(time));
amortisation=schedule(:,3);
prepayment=nan(size(time));
PPC=nan(size(time));
forecast=nan(size(time));
for i=1:size(time,1)
    AO(i)=sum(F(i)*(~tau2(:,i)));
    A(i)=sum(amortisation(i)*(~tau2(:,i)));
    if i>1
        prepayment(i)=AO(i)-AO(i-1)+A(i);
        PPC(i)=-prepayment(i)/AO(i);
        forecast(i)=sum(optimality(~tau2(:,i),i))*(1-exp(-dt*lambda2))*F(i)/AO(i);
    end
end

% Plot
yyaxis left
plot(PPC)
yyaxis right
plot(r)



