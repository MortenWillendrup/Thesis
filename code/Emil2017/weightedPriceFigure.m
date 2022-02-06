%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Price figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
endDate=datenum('2017-06-01');
mkt=market('date',endDate);

% Parameters
F0=100;
n=4;
T=30;
lambda=0.5169;
Alpha=0.1107;
Beta=2.3877;

% Load estimation data
data=load('/users/helmig/dropbox/ku/speciale/Data/estimationData');
info=data.info;
datesXstar=data.dateVector;
XstarCube=data.Xstar;
LAMBDA=data.LAMBDA;

% Load prepaymetns
isin=mkt.bonds.isins{2};
prepay=mkt.bonds.(isin).prepayments(1:end-2);
dates=mkt.bonds.(isin).datesPrepayments(1:end-2);
mat=mkt.bonds.(isin).maturity;

% Load prices
datesPrices=mkt.bonds.(isin).datesPrices;
prices=mkt.bonds.(isin).prices;
B=datesPrices>endDate;
datesPrices(B)=[];prices(B)=[];
time=(datenum(mat)-datesPrices)/365.25;
N=size(time,1);

% Load coupon
R=mkt.bonds.(isin).coupon;

% Pick out bond
B=info(:,2)==mkt.bonds.(isin).coupon&info(:,1)==datenum(mkt.bonds.(isin).maturity);
XstarMat=XstarCube(:,:,B);
Xstar=nan(size(XstarMat,1),1);
for i=1:size(Xstar,1)
    Xstar(i)=hermiteInterpolation(LAMBDA,XstarMat(i,:)',lambda);
end

% Load 6M Cibor
[x,y]=mkt.timeSeries('cibor06M',datesPrices,datesPrices(end));
y=fillmissing(y,'previous');
shortRate=nan(size(datesPrices));
for i=1:size(datesPrices,1)
    B=x==datesPrices(i);
    shortRate(i)=y(B);
end

% Crop Xstar to prepayment dates
B=ismember(datesXstar,dates);
Xstar(~B)=[];

% Initiate model
model=cStantonHW('F0',F0,...
    'R',R,...
    'n',n,...
    'T',T,...
    'lambda',lambda,...
    'A',Alpha,...
    'B',Beta,...
    'market',mkt,...
    'timesteps',3,...
    'spacesteps',400,...
    'rannacher',true,...
    'smoothing',true);


% Set date of market
i=5979;
mkt.date=datesPrices(i);
model.market=mkt;

% Set maturity
model.maturity=floor(time(i)*n)/n;

% Check how far we are
c=sum(dates<=datesPrices(i));

% Get price
model.Xstar=Xstar(1:c);
X=(0:0.01:1)';
price=nan(size(X));
for j=1:size(X,1)
    model.X=X(j);
    [M_a,~,space]=model.pricingFD;
    price(j)=hermiteInterpolation(space,M_a,shortRate(i));
end

% get borrower dist
x=(0:0.001:1)';
cdf=nan(size(x));
pdf=nan(size(x));
for j=1:size(x)
    [cdf(j),pdf(j)]=borrowerDistribution(x(j),Xstar,lambda,Alpha,Beta,4);
end


% Draw a figure
f=figure('color','w','position',[360   305   600   240]);

% Plot prepayments
yyaxis left;
plot(x*100,pdf,...
    'color',[0 0 0])
ylabel('Density function $f_t^X$','interpreter','latex')

% Plot swap
yyaxis right;
p=plot(X*100,price,...
    'color',[0 0 0]);
ylabel('Price','interpreter','latex')
xlabel('Prepayment costs, X (\%)','interpreter','latex')
p(1).LineStyle='--';
p(1).Color=[1 1 1]*0;

% Layout
ax=gca;
for i=1:2
    ax.YAxis(i).TickLabelInterpreter='latex';
    ax.YAxis(i).Label.Interpreter='latex';
    ax.YAxis(i).TickLength=[0,0];
    ax.YAxis(i).Label.FontSize=14;
    ax.YAxis(i).FontSize=14;
    ax.YAxis(i).Color=[0 0 0];
end
ax.XAxis.TickLabelInterpreter='latex';
ax.XAxis.Label.Interpreter='latex';
ax.XAxis.TickLength=[0,0];
ax.XAxis.Label.FontSize=14;
ax.XAxis.FontSize=14;

% Legend stuff
legend({'Implied borrower distribution\quad',...
    'Prices'},...
    'interpreter','latex',...
    'location','northoutside',...
    'orientation','horizontal',...
    'fontsize',12,...
    'box','off')
saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\weightedPriceFigure')