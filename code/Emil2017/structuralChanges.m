%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Structural changes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
mkt=market;
isin='DK0009274300';%isin='DK0009757296';
mat=datenum(mkt.bonds.(isin).maturity);
R=mkt.bonds.(isin).coupon;

% Get housing data
[a,b]=xlsread('/users/helmig/dropbox/ku/speciale/Data/prisudvikling');
houseDates=datenum('1900-01-01')+a(:,1)-1;
houseIndex=[mean(a(:,[2 4]),2),mean(a(:,[3 5]),2)];
houseIndex=houseIndex./repmat(houseIndex(1,:),size(houseIndex,1),1)*100;

% Plot
f=latexPlot('x',houseDates,...
    'y',houseIndex,...
    'legend',{'Single Family Houses\quad\quad';'Apartments'},...
    'location','northoutside',...
    'ylabel','Index',...
    'orientation','horizontal');datetick;
f.Children(2).XTick=datenum(strcat(cellstr(num2str((1993:2:2017)')),'-01-01'));
f.Children(2).XTickLabel=cellstr(num2str((1993:2:2017)'));
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\houseindex')

% Bond prices
bondDates=(mkt.bonds.(isin).datesPrices(1):mkt.bonds.(isin).datesPrices(end))';
bondPrices=nan(size(bondDates));
for i=1:size(bondPrices,1)
    B=mkt.bonds.(isin).datesPrices==bondDates(i);
    if any(B)
        bondPrices(i)=mkt.bonds.(isin).prices(B);
    end
end
bondPrices=fillmissing(bondPrices,'previous');

% get schedule
startDate=floor(mat-31.5*365.25);
scedule=annuity(R,80,30,4);
temp=startDate+floor(scedule(:,1)*365.25)+10;
scedule(:,1)=datenum(strcat(num2str(year(temp)),'-',num2str(month(temp)),'-01'));

% Remove prices before start dates
B=bondDates<startDate;
bondDates(B)=[];bondPrices(B)=[];

% Normalise bond prices
P=min(bondPrices./100,1);

% Loan value
loanNominal=nan(size(P));
for i=1:size(P,1)
    B=scedule(:,1)==bondDates(i);
    if any(B)
        loanNominal(i)=scedule(B,5);
    end
end
loanNominal=fillmissing(loanNominal,'previous');
loanNominal(isnan(loanNominal))=80;


% House value
houseVal=nan(size(P));
for i=1:size(houseVal,1)
    B=houseDates<=bondDates(i);
    houseVal(i)=houseIndex(sum(B),1);
end
houseVal=houseVal./houseVal(1);

temp=loanNominal./houseVal;
% B=[true;~(diff(temp)==0)];
% plot(bondDates(B),temp(B))
% ylabel('LTV')
% yyaxis right
% plot(bondDates(1:7:end),bondPrices(1:7:end))
% ylabel('Price')
% datetick


if true
    % Parameters
    F0=100;
    n=4;
    T=30;
    lambda=0.5773;
    Alpha=0.0612;
    Beta=2.2511;

    % Load estimation data
    data=load('/users/helmig/dropbox/ku/speciale/Data/estimationData');
    info=data.info;
    datesXstar=data.dateVector;
    XstarCube=data.Xstar;
    LAMBDA=data.LAMBDA;

    % Load prepaymetns dates
    dates=mkt.bonds.(isin).datesPrepayments(1:end-2);
    mat=mkt.bonds.(isin).maturity;

    % Load prices
    endDate=datenum('2013-12-31');
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
        'timesteps',1,...
        'spacesteps',300,...
        'rannacher',false,...
        'smoothing',false);

    % Loop over dates and get theoretical prices
    pricesTheo=nan(size(prices));c=0;
    for i=1:7:size(datesPrices,1)
        % Set date of market
        mkt.date=datesPrices(i);
        model.market=mkt;

        % Set maturity
        model.maturity=floor(time(i)*n)/n;

        % Check how far we are
        c=sum(dates<=datesPrices(i));
        
        % Get price
        model.Xstar=Xstar(1:c);
        pricesTheo(i)=model.getPrice(shortRate(i));%,prices(i));
    end

end



% Draw a figure
f=figure('color','w','position',[360   305   600   273]);

% Plot
B=[true;~(diff(temp)==0)];
p=plot(bondDates(B),temp(B).*P(B)./P(1),...
    'color',[0 0 0]);
ylabel('LTV (\%)','interpreter','latex')
xlabel('Year','interpreter','latex')
p(1).LineStyle='-';
p(1).Color=[0 0 0];
p(1).LineWidth=2;
ylim([75 105])

yyaxis right
plot(bondDates(1:7:end),bondPrices(1:7:end),...
    'color',[0 0 0]*0.65,...
    'linewidth',1)
hold on
plot(datesPrices(1:7:end),pricesTheo(1:7:end),...
    'color',[1 1 1]*0.65,...
    'linestyle','-',...
    'linewidth',2)

ylabel('Price','interpreter','latex')
datetick

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
legend({'LTV (Left)\quad',...
    'Market Price (Right)\quad',...
    'Model Price (Right)\quad'},...
    'interpreter','latex',...
    'location','northoutside',...
    'orientation','horizontal',...
    'fontsize',10.8,...
    'box','off')
xlim([datenum('2006-06-30') datenum('2013-06-01')])
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\LTV')


isins=mkt.bonds.isins;
s=nan(10,1);
for i=48:58
    s(i-47)=max(mkt.bonds.(isins{i}).outstanding);
end
sum(s)/1e9