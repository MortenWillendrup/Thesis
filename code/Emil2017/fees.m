%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fees
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Market
mkt=market('date','2017-06-01');

% Model Setting
Alpha=0.0612;
Beta=2.2511;
lambda=0.5773;
R=0.025;
mat=30;
F0=100;
feeRateLong=0.008;
feeRateShort=0.008;
model=cStantonHW('market',mkt,...
    'A',Alpha,...
    'B',Beta,...
    'lambda',lambda,...
    'maturity',mat,...
    'coupon',R,...
    'notional',F0,...
    'FeeRateLong',feeRateLong,...
    'FeeRateShort',feeRateShort,...
    'PrepaymentCosts',0.15,...
    'terms',4);

% Get values
[Assets,Liabilities,space]=model.pricingFD;

% Adjust fees
model.feeRateLong=0.0065;
model.feeRateShort=0.095;

% New valus
[Assets2,Liabilities2]=model.pricingFD;

% Plot
f=figure('color','w','position',[360   494   594   220]);
subplot(1,2,1)
latexPlot('x',space*100,...
    'y',[Assets,Assets2],...
    'f',gca,...
    'xlabel','Short rate (\%)',...
    'ylabel','Price',...
    'legend',{'Before';'After'})
xlim([-2,5])
ax=gca;
ax.Position(2)=ax.Position(2)+0.08;
ax.Position(4)=ax.Position(4)-0.08;
subplot(1,2,2)
latexPlot('x',space*100,...
    'y',Assets2-Assets,...
    'f',gca,...
    'xlabel','Short rate (\%)',...
    'ylabel','Price change')
ax=gca;
ax.Position(2)=ax.Position(2)+0.08;
ax.Position(4)=ax.Position(4)-0.08;
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\assetChange')


% Fee rates over time
[a,b,c]=xlsread('/users/helmig/dropbox/ku/speciale/Data/bidragssatser.xlsx');
feeRates=a(:,[2 3 6]);
names=c(1,[2 3 6])';
dates=a(:,1)+datenum('1900-01-01');
B=dates<=datenum('2008-01-01');
dates(B)=[];
feeRates(B,:)=[];

l={'\quad Average level\quad';'\quad$\leq$ 1 Year\quad';'\quad$\geq$ 10 Year'};
% Figure
f=latexPlot('x',dates,...
    'y',feeRates,...
    'legend',l,...
    'location','northoutside',...
    'orientation','horizontical',...
    'ylabel','Fee rate (\%)');
datetick;
ylim([0.2 1.2])
grid on
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\feeRates')

