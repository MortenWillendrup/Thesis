%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Figures
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
addpath('/Users/helmig/Dropbox/KU/Speciale/Data')

% Load market
mkt=market('date','2017-06-01');

% Extract rates
startDate='2000-01-01';
endDate='2017-06-01';
[dateVector,rates]=timeSeries(mkt,'CIBOR01M',startDate,endDate);
rates=fillmissing(rates,'previous');

% Plot
f=latexPlot('dates',dateVector,...
          'y',rates*100,...
          'ylabel','',...
          'position',[360 278 600 300]);

% Fix dates for x ticks
dateMatrix=datevec(dateVector);
B=dateMatrix(:,2)==1&dateMatrix(:,3)==1;
xDateTicks=dateVector(B);
xDateTicks=xDateTicks(1:4:end);
xDateString=num2str(datestr(xDateTicks,10));
ax=gca;
ax.XAxis.TickValues=xDateTicks;
ax.XAxis.TickLabels=xDateString;
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\swaprates')

%% Swap calibration
close all
B=mkt.ciborTenors<=0.5;

% Calibrate market data
[zero,forward,zeroTenors]=swapCalibration(mkt.swapTenors,...
                    mkt.swapCurve,...
                    mkt.ciborTenors(B),...
                    mkt.ciborCurve(B));

% Calculate interpolated zero rates
tenors=(min(mkt.ciborTenors):0.01:max(mkt.swapTenors))';
[zeroRates,fwdRates]=hermiteInterpolationFwd(zeroTenors,zero,tenors);

f=latexPlot('x',tenors,...
            'y',[zeroRates,fwdRates]*100,...
            'xlabel','Maturity',...
            'ylabel','$\%$');
hold on
scatter(mkt.swapTenors,mkt.swapCurve*100,...
          'marker','o',...
          'MarkerEdgeColor',[0 0 0]);
scatter(mkt.ciborTenors(B),mkt.ciborCurve(B)*100,...
          'marker','o',...
          'MarkerEdgeColor',[0 0 0]);

legendStr={'Zero Curve';'Forward Curve';'Cibor/Swap rates'};
legend(legendStr,'fontsize',11,...
    'interpreter','latex',...
    'EdgeColor',[1,1,1],...
    'location','northoutside',...
    'Orientation','horizontal');
set(f,'position',[360   221   600   307])
xlim([-0.5 30.5])
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\swapCurve')

%% Cap calibration
close all

% Set date
mkt.date=datenum('2017-06-01');

% Calibrate swap
B=mkt.ciborTenors<=0.5;
[zero,~,zeroTenors]=swapCalibration2(mkt.swapTenors,...
                    mkt.swapCurve/100,...
                    mkt.ciborTenors(B),...
                    mkt.ciborCurve(B)/100);
                
% Calibrate CAP
[kappa,sigma]=capCalibration(mkt.capTenors,mkt.capCurve,zeroTenors,zero);

% Plot prices against theoretical prices
capTenorsTheo=(2:21)';
capPrices=ATMcapPriceHW(capTenorsTheo,kappa,sigma,zeroTenors,zero)*10000;

f=latexPlot('x',capTenorsTheo,...
            'y',capPrices,...
            'xlabel','Maturity',...
            'ylabel','BPS',...
            'position',[360   305   600   273]);
hold on
scatter(mkt.capTenors,mkt.capCurve,...
          'marker','o',...
          'MarkerEdgeColor',[0 0 0]);
xlim([0 22])
legend({'Model Quotes';'Market Quotes'},...
    'fontsize',11,...
    'interpreter','latex',...
    'EdgeColor',[1,1,1],...
    'location','northoutside',...
    'Orientation','horizontal');
saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\capCalibration')

%% Kappa and Sigma over time
dateVector=(datenum('2002-01-01'):7:datenum('2017-06-01'))';
SIGMA=nan(size(dateVector));
h=waitbar(0);
N=size(dateVector,1);
for i=1:N
    try
        % Waitbar
        waitbar(i/N,h)
        
        % Date
        mkt.date=dateVector(i);
                               
        % CAP calibration        
        SIGMA(i)=capCalibrationKnownKappa(mkt,kappa);
    catch 
    end
   
end

%% Plot
ydata=fillmissing(SIGMA,'previous')*100;
f=latexPlot('x',dateVector,...
            'y',ydata,...
            'ylabel','Volatility \%');

% Fix dates for x ticks
dateVector2=(dateVector(1):dateVector(end))';
dateMatrix=datevec(dateVector2);
B=dateMatrix(:,2)==1&dateMatrix(:,3)==1;
xDateTicks=dateVector2(B);
xDateTicks=xDateTicks(1:2:end);
xDateString=num2str(datestr(xDateTicks,10));
ax=gca;
ax.XAxis.TickValues=xDateTicks;
ax.XAxis.TickLabels=xDateString;
xlim([dateVector2(1)+50,dateVector2(end)])
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\kappaSigma')


%% PRepayment figure
prepayments=mkt.bonds.DK0009753469.Extraordinary;
dates=mkt.bonds.DK0009753469.Dates;
[swapDates,swap]=mkt.timeSeries('swap10y',dates(1),dates(end));
swap=fillmissing(swap,'previous');

% Plot prepayments
f=figure('color','w','position',[360   305   600   273]);
yyaxis left;
b=bar(dates,prepayments*100,...
      'barWidth',1,...
      'FaceColor',[0.8,0.8,0.8]);
ylabel('Prepayments (\%)','interpreter','latex')

% Plot costs
yyaxis right;
plot(swapDates(1:7:end),swap(1:7:end)*100,...
       'color',[0 0 0])
ylabel('10Y Swap Rate (\%)','interpreter','latex')
xlabel('Years','interpreter','latex')

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
legend({'Prepayments (Left)\quad',...
        '10Y Swap Rate (Right)'},...
        'interpreter','latex',...
        'location','northoutside',...
        'orientation','horizontal',...
        'fontsize',12,...
        'box','off')
    
% Fix dates for x ticks
dateVector=(swapDates(1):swapDates(end))';
dateMatrix=datevec(dateVector);
B=dateMatrix(:,2)==1&dateMatrix(:,3)==1;
xDateTicks=dateVector(B);
xDateTicks=xDateTicks(2:4:end);
xDateString=num2str(datestr(xDateTicks,10));
ax.XAxis.TickValues=xDateTicks;
ax.XAxis.TickLabels=xDateString;
xlim([swapDates(1)-90 swapDates(end)+90])


% Save
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\DK0009753469')

