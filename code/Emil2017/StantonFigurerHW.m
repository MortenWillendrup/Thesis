%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Extended Stanton - figurer
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc;clear;close all;

% Set parameters
kappa=0.0796;
X=0.15;
F0=100;
R=0.035;
n=4;
T=30;
lambda=1.5;

% Load market
mkt=market('date','2017-06-01');

% Initiate model
model=cStantonHW('kappa',kappa,...
                 'X',X,...
                 'F0',F0,...
                 'R',R,...
                 'n',n,...
                 'T',T,...
                 'lambda',lambda,...
                 'market',mkt);

%% Loop over different rho (Stanton fig 5)
if false
    
tic
%lambda=[0,0.1,0.3,0.6,2,10000];
lambda=[0.3,0.6,0.9,1.2,1.5,1.8];
prices=[];
for i=1:max(size(lambda))
    model.lambda=lambda(i);
    [M_a,~,space]=model.pricingFD;
    prices=[prices,M_a];
end
toc

% Plot prices from r=0% to r=25%
b=space>=-0.05&space<=10;
f=latexPlot('x',space(b)*100,...
          'y',prices(b,:),...
          'xlabel','Short Rate (\%)',...
          'ylabel','Price',...
          'legend',cellstr(strcat('$\lambda_2=',num2str(lambda'),'$')),...
          'position',[360 278 600 300]);
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\FD1HW')
   

%% Loop over different rho and X(Stanton fig 5)
tic
lambda=[0.75 1.5];
X=[0.1 0.2];
prices=[];
leg={};
for j=1:2
    for i=1:2
        model.lambda=lambda(i);
        model.X=X(j);
        [M_a,~,space]=model.pricingFD;
        prices=[prices,M_a];
        leg=[leg;{sprintf('$\\lambda=%s, X=%s$',...
                  num2str(lambda(i),'%4.2f'),num2str(X(j),'%4.2f'))}];
    end
end
toc

% Plot prices from r=0% to r=25%
b=space>=-0.055&space<=0.075;
f=latexPlot('x',space(b)*100,...
          'y',prices(b,:),...
          'xlabel','Short Rate (\%)',...
          'ylabel','Price',...
          'legend',leg,...
          'position',[360 278 600 240]);
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\FD1HW')
   

%% Loop over different coupons
coupon=[0.02,0.025,0.03,0.04,0.05,0.06];
prices=[];
model.lambda=1.5;
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
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\FD2HW')
   
%% Loop over different X (Stanton fig 6)
X=[0,7,16,24,40,100]/100;
prices=[];
model.lambda=1.5;
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
      
%% Delivery option or not
model.lambda=2.5;
model.X=0.20;
model.delivery=true;
[M_a,M_l]=model.pricingFD;
model.delivery=false;
[M_a_no_delivery,M_l_no_delivery,space]=model.pricingFD;
assets=[M_a,M_a_no_delivery];
liabilities=[M_l,M_l_no_delivery];

% Plot prices from r=0% to r=25%
b=space<=0.25;
f=latexPlot('x',space(b)*100,...
          'y',[assets(b,:),liabilities(b,:)],...
          'xlabel','Short Rate (\%)',...
          'ylabel','Price',...
          'legend',{'Assets';'Assets No del.';'Liabilities';'Liabilities No del.'},...
          'position',[360 278 600 300]);
      
end
%% Plot against swap rate
close all
model.lambda=2;
model.X=0.1;
[dateVector,swap10y]=mkt.timeSeries('swap10y','2002-01-01','2017-06-01');
[~,cibor1M]=mkt.timeSeries('cibor01m','2002-01-01','2017-06-01');
[~,cibor6M]=mkt.timeSeries('cibor06m','2002-01-01','2017-06-01');
swap10y=fillmissing(swap10y,'previous');
cibor1M=fillmissing(cibor1M,'previous');

% Pre allocation
bondPrices=nan(size(dateVector,1),1);
impliedVol=nan(size(dateVector,1),1);

% Waitbar
h=waitbar(0);

% Loop to generate data
for i=1:30:size(dateVector,1)
    mkt.date=dateVector(i);
    model.market=mkt;
    [Ma,~,r]=model.pricingFD;
    bondPrices(i)=hermiteInterpolation(r,Ma,cibor1M(i));
    impliedVol(i)=model.sigma;
    waitbar(i/size(dateVector,1),h)
end
close(h)

f=figure();
B=~isnan(bondPrices);
yyaxis left;plot(dateVector(B),bondPrices(B))
yyaxis right;plot(dateVector(B),swap10y(B))


