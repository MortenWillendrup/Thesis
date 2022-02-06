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
lambda=0.5773;
Alpha=0.0612;
Beta=2.2511;

% Gamma numbers
% lambda=0.1658;
% Alpha=0.2230;
% Beta=0.6000;
% 
% 
% Alpha=-4.6050;
% Beta=1.8896;
% lambda=0.5890;

% Load estimation data
data=load('/users/helmig/dropbox/ku/speciale/Data/estimationData');
info=data.info;
datesXstar=data.dateVector;
XstarCube=data.Xstar;
LAMBDA=data.LAMBDA;

% Load prepaymetns
isin=mkt.bonds.isins{2};isin='DK0009761645';isin='DK0009757296';
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

% Find prices above par
abovePar=true(size(Xstar,1),1);
for i=1:size(Xstar,1)
    B=datesPrices<=datesXstar(i);
    if any(B)
        if prices(sum(B))<100
            abovePar(i)=false;
        end
    end
end

% Crop Xstar to prepayment dates
B=ismember(datesXstar,dates);
Xstar(~B)=[];
abovePar(~B)=[];

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

% Loop over dates and get theoretical prices
pricesTheo=nan(size(prices));
OAS=nan(size(prices));
OAD=nan(size(prices,1),2);
OAC=nan(size(prices,1),2);
Carry=nan(size(prices,1),2);
c=0;
for i=1:100:size(datesPrices,1)
    % Set date of market
    mkt.date=datesPrices(i);
    model.market=mkt;

    % Set maturity
    model.maturity=floor(time(i)*n)/n;

    % Check how far we are
    c=sum(dates<=datesPrices(i));

    % Get price
    model.Xstar=Xstar(1:c);
    pricesTheo(i)=model.getPrice(shortRate(i));
    [OAD(i,1),OAC(i,1),~,~,Carry(i,1)]=model.keyfigures(shortRate(i));
    %OAS(i)=model.getOAS(shortRate(i),prices(i));
    %[OAD(i,2),OAC(i,2),~,~,Carry(i,2)]=model.keyfigures(shortRate(i));
end

% Risk source
s=load('/users/helmig/dropbox/ku/speciale/Data/yields','dateVector','yields','forwards');
dateVector=s.dateVector;
yields=s.yields;
forwards=s.forwards;
dr=diff(shortRate(1:7:end));
y5=nan(size(datesPrices));
y10=nan(size(datesPrices));
f10=nan(size(datesPrices));
y15=nan(size(datesPrices));
for i=1:size(datesPrices,1)
    B=dateVector==datesPrices(i);
    try
        y5(i)=yields(B,4);
        y10(i)=yields(B,5);
        y15(i)=yields(B,6);
    catch
    end
end
y5=fillmissing(y5,'previous');
y10=fillmissing(y10,'previous');
y15=fillmissing(y15,'previous');
B=1/mkt.kappa*(1-exp(-mkt.kappa*10));
dr3=10/B*(diff(y10(1:7:end))+(y15(8:7:end)-y5(8:7:end))./10*1/52);

% Test of keyfigures
Duration=OAD(1:7:end-7,2);
Conv=OAC(1:7:end-7,2);
TimeValue=Carry(1:7:end-7,2);
hedgeSeries=prices(1)+[0;cumsum(TimeValue*1/52 ...
    -Duration.*dr*100 ...
    +1/2*Conv.*dr.^2*10000)];

% Test of keyfigures
Duration=OAD(1:7:end-7,2);
Conv=OAC(1:7:end-7,2);
TimeValue=Carry(1:7:end-7,2);
hedgeSeries2=prices(1)+[0;cumsum(TimeValue*1/52 ...
    -Duration.*dr3*100 ...
    +1/2*Conv.*dr3.^2*10000)];

% Draw a figure
f=figure('color','w','position',[360   305   600   273]);

% Prices
datesToPlot=datesPrices(1:7:end);
pricesToPlot=[prices(1:7:end),hedgeSeries2];
p=plot(datesToPlot,pricesToPlot,...
    'color',[0 0 0]);
ylabel('Price','interpreter','latex')
xlabel('Year','interpreter','latex')
p(1).LineStyle='-';
p(2).LineStyle='-';
p(1).Color=[1 1 1]*0.65;
p(1).LineWidth=2;

% Layout
ax=gca;
for i=1:1
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
legend({'Market Price\quad',...
    'Hedge Series (using $y_t^{t+10}$)\quad'},...
    'interpreter','latex',...
    'location','northoutside',...
    'orientation','horizontal',...
    'fontsize',10.8,...
    'box','off')

% Fix dates for x ticks
dateVector=(datesToPlot(1):datesToPlot(end))';
dateMatrix=datevec(dateVector);
B=dateMatrix(:,2)==1&dateMatrix(:,3)==1;
xDateTicks=dateVector(B);
xDateTicks=xDateTicks(1:2:end);
xDateString=num2str(datestr(xDateTicks,10));
ax.XAxis.TickValues=xDateTicks;
ax.XAxis.TickLabels=xDateString;
xlim([datesToPlot(1)-90 datesToPlot(end)+90])

%saveFigAsPdf(f,strcat('\users\helmig\dropbox\KU\speciale\thesis\hedge10Y',isin))

% mats=(30:-0.25:0.25)';
% model.Xstar=[];
% priser=nan(size(mats));
% OAD2=nan(size(mats));
% OAC2=nan(size(mats));
% Carry2=nan(size(mats));
% h=waitbar(0);
% for i=1:size(mats,1)
%     waitbar(i/size(mats,1),h)
%     % Set maturity
%     model.maturity=mats(i);
%     priser(i)=model.getPrice(0);
%     [OAD2(i,1),OAC2(i,1),~,~,Carry2(i,1)]=model.keyfigures(0);
%     fprintf('\n%4.2f',Carry2(i,1))
% end
% close(h)

% Draw a figure
f=figure('color','w','position',[360   305   600   273]);

% Plot prepayments
yyaxis left;
b=bar(dates,prepay,...
    'barWidth',1,...
    'FaceColor',[0.8,0.8,0.8],...
    'facealpha',1);
hold on
modelPrepay=modelPrepayments(Alpha,Beta,lambda,Xstar,n)*100;
modelPrepay(~abovePar)=0;
plot(dates,modelPrepay,...
    'color',[0 0 0],...
    'linewidth',2)
ylabel('Prepayments','interpreter','latex')

% Plot swap
yyaxis right;
datesToPlot=datesPrices(1:7:end);
B=(datesToPlot<datenum('1993-12-31'));datesToPlot(B)=[];
pricesToPlot=[prices(1:7:end),pricesTheo(1:7:end)];pricesToPlot(B,:)=[];
p=plot(datesToPlot,pricesToPlot,...
    'color',[0 0 0]);
ylabel('Price','interpreter','latex')
xlabel('Year','interpreter','latex')
p(1).LineStyle='-';
p(2).LineStyle='-';
p(1).Color=[1 1 1]*0.65;
p(1).LineWidth=2;

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
legend({'Prepayments\quad',...
    'Estimated Prepayments\quad',...
    'Market Price\quad',...
    'Model Price'},...
    'interpreter','latex',...
    'location','northoutside',...
    'orientation','horizontal',...
    'fontsize',10.8,...
    'box','off')

% Fix dates for x ticks
dateVector=(datesToPlot(1):datesToPlot(end))';
dateMatrix=datevec(dateVector);
B=dateMatrix(:,2)==1&dateMatrix(:,3)==1;
xDateTicks=dateVector(B);
xDateTicks=xDateTicks(1:2:end);
xDateString=num2str(datestr(xDateTicks,10));
ax.XAxis.TickValues=xDateTicks;
ax.XAxis.TickLabels=xDateString;
xlim([datesToPlot(1)-90 datesToPlot(end)+90])

%saveFigAsPdf(f,strcat('\users\helmig\dropbox\KU\speciale\thesis\price',isin))


% Comets
% f=figure;hold on;
% xlim([0.00 0.09]);ylim([60 130]);
% for i=1:7:size(swap10y,1)
%     s=scatter(swap10y(i),pricesTheo(i));s.MarkerEdgeColor=[0.2 0.2 0.8];
%     s=scatter(swap10y(i),prices(i));s.MarkerEdgeColor=[0.6 0.3 0.1];
%
%     pause(0.01)
% end