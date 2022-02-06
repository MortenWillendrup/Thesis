%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   price scatter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Market
mkt=market('date','2017-06-01');
isin=mkt.bonds.isins{74}; % NYK 4'41 DK0009775355 - 74

% Get bond prices
prices=mkt.bonds.(isin).prices;
priceDates=mkt.bonds.(isin).datesPrices;

% Get swap
[dates,swap]=mkt.timeSeries('swap10y',priceDates(1),priceDates(end));

% Match
temp=nan(size(swap));
for i=1:size(swap,1)
    B=priceDates==dates(i);
    if any(B)
        temp(i)=prices(B);
    end
end
prices=temp;


% Draw a figure
f=figure('color','w','position',[360   305   600   273]);

% Scatter
yearsToScatter=[2010 2011;
                2012 2013;
                2014 2015;
                2016 2017];
yearVector=year(dates);
for i=1:size(yearsToScatter,1)
    B=ismember(yearVector,yearsToScatter(i,:)');
    s=scatter(swap(B)*100,prices(B),'markerfacecolor',[1 1 1]*i/5,...
        'markeredgecolor',[0 0 0]);
    s.MarkerEdgeAlpha=0.2;
    hold on
end
ylabel('Price','interpreter','latex')
xlabel('10Y Swap (\%)','interpreter','latex')


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
legend(strcat('\textit{',cellstr(num2str(yearsToScatter(:,1))),'-',cellstr(num2str(yearsToScatter(:,2))),'}'),...
    'interpreter','latex',...
    'location','northoutside',...
    'orientation','horizontal',...
    'fontsize',10.8,...
    'box','off')
box on
saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\scatterPrices')
