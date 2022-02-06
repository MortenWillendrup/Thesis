%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   First Pricing example
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;

% Set parameters
kappa=0.29368;
mu=0.07935;
sig=0.11425;
q=-0.12165;
X=0.24;
F0=100;
R=0.125;
n=12;
T=30;
lambda1=0;%0.0338;
lambda2=0.6452;

% Initiate model
model=cStanton('kappa',kappa,...
               'mu',mu,...
               'sigma',sig,...
               'q',q,...
               'X',X,...
               'F0',F0,...
               'R',R,...
               'n',n,...
               'T',T,...
               'lambda1',lambda1,...
               'lambda2',lambda2);

%% Loop over different rho (Stanton fig 5)
tic
RHO=[0,0.1,0.3,0.6,2,10000];
prices=[];
for i=1:6
    model.lambda2=RHO(i);
    [M_a,~,space]=model.pricingFD;
    prices=[prices,M_a];
end
toc
% Plot prices from r=0% to r=25%
b=space<=0.25;
f=latexPlot('x',space(b)*100,...
          'y',prices(b,:),...
          'xlabel','Short Rate (\%)',...
          'ylabel','Price',...
          'legend',cellstr(strcat('$\lambda_2=',num2str(RHO'),'$')),...
          'position',[360 278 600 300]);
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\FD1')
   

%% Loop over different coupons
coupon=[0.02,0.03,0.05,0.075,0.1, 0.125];
prices=[];
model.lambda2=2;
for i=1:6
    model.R=coupon(i);
    [M_a,~,space]=model.pricingFD;
    prices=[prices,M_a];
end

% Plot prices from r=0% to r=25%
b=space<=0.25;
f=latexPlot('x',space(b)*100,...
          'y',prices(b,:),...
          'xlabel','Short Rate (\%)',...
          'ylabel','Price',...
          'legend',cellstr(strcat('$R=',num2str(100*coupon'),'\%$')),...
          'position',[360 278 600 300]);
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\FD2')
   
%% Loop over different X (Stanton fig 6)
X=[0,7,16,24,40,100]/100;
prices=[];
model.lambda1=0;
model.lambda2=1000000;
for i=1:6
    model.X=X(i);
    [M_a,~,space]=model.pricingFD;
    prices=[prices,M_a];
end

% Plot prices from r=0% to r=25%
b=space<=0.25;
f=latexPlot('x',space(b)*100,...
          'y',prices(b,:),...
          'xlabel','Short Rate (\%)',...
          'ylabel','Price',...
          'legend',cellstr(strcat('$X=',num2str(X'*100),'\%$')),...
          'position',[360 278 600 300]);
%% Loop over different X with estimated parameters (Stanton fig 7)
X=[0,40,100000]/100;
prices=[];
model.lambda1=0.0338;
model.lambda2=0.6452;
for i=1:3
    model.X=X(i);
    [M_a,~,space]=model.pricingFD;
    prices=[prices,M_a];
end

% Plot prices from r=0% to r=25%
b=space<=0.25;
plot(space(b),prices(b,:))
legend(cellstr(strcat('$X=',num2str(X'*100),'\%$')),...
       'Interpreter','latex',...
       'FontSize',12)

%% Check againts CIR closed form
model.lambda1=0;
model.lambda2=0;
[M_a,~,space]=model.pricingFD;

schedule=annuity(R,F0,T,n);
k=kappa+q;
theta=kappa*mu/k;
h=sqrt(k^2+2*sig^2);
A=@(T)(2*h*exp((h+k)*T/2)/(2*h+(h+k)*(exp(h*T)-1)))^(2*k*theta/sig^2);
B=@(T)2*(exp(h*T)-1)/(2*h+(h+k)*(exp(h*T)-1));

truePrices=nan(size(space));
for j=1:size(space,1)
    P=nan(size(schedule,1),1);
    for i=1:size(P,1)
        t=schedule(i,1);
        P(i)=A(t)*exp(-space(j)*B(t));
    end
    truePrices(j)=sum(schedule(:,4).*P);
end
b=space<=0.25;
yyaxis left
p=plot(space(b),truePrices(b),'LineWidth',2);
ylabel('Price')
hold on
p=scatter(space(b),M_a(b));
yyaxis right
plot(space(b),M_a(b)-truePrices(b))
legend({'True price';'Finite Difference';'Deviation';'123'})
xlabel('$r$','Interpreter' ,'latex')
ylabel('Deviation')
set(gca,'fontsize',14)



%{ 
%%% Check againts CIR closed form
hold on
k=kappa+q;
theta=kappa*mu/k;
h=sqrt(k^2+2*sig^2);
A=@(T)(2*h*exp((h+k)*T/2)/(2*h+(h+k)*(exp(h*T)-1)))^(2*k*theta/sig^2);
B=@(T)2*(exp(h*T)-1)/(2*h+(h+k)*(exp(h*T)-1));
P=nan(size(space));
for i=1:size(space,1)
    P(i)=A(T)*exp(-space(i)*B(T));
end
plot(space,P)
%}


