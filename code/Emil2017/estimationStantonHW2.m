%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Stanton estimation under Hull & White dynamics (multiple bonds)
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
endDate=datenum('2017-06-01');
mkt=market('date',endDate);

% Estimation date
estimationDate=datenum('2017-06-01');

% Load estimation data
data=load('/users/helmig/dropbox/ku/speciale/Data/estimationData');
dateVector=data.dateVector;
info=data.info;
Xstar=data.Xstar;
LAMBDA=data.LAMBDA;

% Load prepaymetns
isins=mkt.bonds.isins;
N=size(dateVector,1);
N2=size(LAMBDA,1);
N3=size(isins,1);
prepayments=nan(N,N3);
maturities=nan(N3,1);
coupons=nan(N3,1);
time=nan(N,N3);
XstarMatrix=nan(N,N2,N3);
abovePar=ones(N,N3);
for i=1:N3
    % Get maturity 
    coupons(i)=mkt.bonds.(isins{i}).coupon;
    maturities(i)=datenum(mkt.bonds.(isins{i}).maturity);
    time(:,i)=round((maturities(i)-dateVector)/365*4,0)/4;
    
    % Load prices
    prices=mkt.bonds.(isins{i}).prices;
    pricesDates=mkt.bonds.(isins{i}).datesPrices;
    for j=1:N
        B=pricesDates<=dateVector(j);
        if any(B)
            if prices(sum(B))<100
                abovePar(j,i)=0;
            end
        end
    end
    
    % Extract prepayments
    temp1=mkt.bonds.(isins{i}).prepayments/100;
    temp2=mkt.bonds.(isins{i}).datesPrepayments;
    %bar(temp1)
    % Remove issuance period
    B=false(size(temp1));
    for j=1:size(temp1,1)
        if temp1(j)<0.0005
            B(j)=true;
        else
            break
        end
    end
    temp1(B)=[];temp2(B)=[];
    
    % Remove dates after estimationDate
    B=temp2>estimationDate;
    temp1(B)=[];temp2(B)=[];
    
    % Correct wrong dates
    temp3=datevec(temp2);
    temp2=temp2+ones(size(temp2))-temp3(:,3);
    
    % Save remaining data
    B=ismember(dateVector,temp2);
    prepayments(B,i)=temp1;
    
    % Extract Xstar and insert into Xstar Matrix
    B=maturities(i)==info(:,1)&coupons(i)==info(:,2);
    XstarMatrix(:,:,i)=Xstar(:,:,B);
end

% Remove unwanted isins
IX=[36 85:116]';
prepayments(:,IX)=[];
XstarMatrix(:,:,IX)=[];
isins(IX)=[];
abovePar(:,IX)=[];

% Load 1M Cibor
[x,y]=mkt.timeSeries('cibor01m',dateVector(1)-10,dateVector(end)+10);
y=fillmissing(y,'previous');
cibor=nan(N,1);
for i=1:N
    B=x==dateVector(i);
    cibor(i)=y(B);
end

% Load 10Y Swap
[x,y]=mkt.timeSeries('swap10Y',dateVector(1)-10,dateVector(end)+10);
y=fillmissing(y,'previous');
swap10y=nan(N,1);
for i=1:N
    B=x==dateVector(i);
    swap10y(i)=y(B);
end

% Estimate
Alpha=3;
Beta=3;
lambda=1;
theta=[Alpha;Beta;lambda];
lb=zeros(3,1);
ub=[100;100;3];
tic
theta=fmincon(@(theta)GMMform2(theta,prepayments,XstarMatrix,LAMBDA,abovePar),theta,[],[],[],[],lb,ub)
toc


% Find prediction errors over time
expPrepay1=nan(size(dateVector,1),size(isins,1));
expPrepay2=nan(size(dateVector,1),size(isins,1));
expPrepayProb1=nan(size(dateVector,1),size(isins,1));
expPrepayProb2=nan(size(dateVector,1),size(isins,1));

% Get expected prepayments
%theta=[0.1776;4.3242;0.6997];
for k=1:size(isins,1)
    B=~isnan(prepayments(:,k));
    % Interpolate Xstar given lambda
    XstarInt=nan(size(XstarMatrix,1),1);
    for i=1:size(Xstar,1)
        XstarInt(i)=hermiteInterpolation(LAMBDA,XstarMatrix(i,:,k)',theta(3));
    end
    [p1,p2]=modelPrepayments(theta(1),theta(2),theta(3),XstarInt(B),4);            
    p1(~abovePar(B,k))=0;
    p2(~abovePar(B,k))=0;
    expPrepay1(B,k)=p1';
    expPrepayProb1(B,k)=p2';
end

% Get expected prepayments
for k=1:size(isins,1)
    B=~isnan(prepayments(:,k));
    % Interpolate Xstar given lambda
    XstarInt=nan(size(XstarMatrix,1),1);
    for i=1:size(Xstar,1)
        XstarInt(i)=hermiteInterpolation(LAMBDA,XstarMatrix(i,:,k)',theta(3));
    end
    [p1,p2]=modelPrepayments(theta(1),theta(2),theta(3),XstarInt(B),4);            
    p1(~abovePar(B,k))=0;
    p2(~abovePar(B,k))=0;
    expPrepay2(B,k)=p1';
    expPrepayProb2(B,k)=p2';
end

% Calc market probs
marketProbs=nan(size(prepayments));
for k=1:size(isins,1)
    temp=prepayments(:,k);
    B=~isnan(temp);
    temp2=temp(B);
    marketRemain=cumprod(1-[0;temp2(1:end-1)]);
    marketProbs(B,k)=temp2.*marketRemain;
end

% Forecast plot
avgPrepay=nan(size(prepayments,1),1);
avgExp1=nan(size(prepayments,1),1);
avgExp2=nan(size(prepayments,1),1);
for i=1:size(avgPrepay,1)
    B=~isnan(prepayments(i,:));
    avgPrepay(i)=1/sum(B)*sum(prepayments(i,B));
    avgExp1(i)=1/sum(B)*sum(expPrepay1(i,B));
    avgExp2(i)=1/sum(B)*sum(expPrepay2(i,B));
end
B1=dateVector<datenum('2011-01-01');
B2=dateVector>=datenum('2010-09-30');

% Plot prepayments
f=figure('color','w','position',[360   305   600   273]);

% Plot prepayments 
b=bar(dateVector,avgPrepay*100,...
      'barWidth',1,...
      'FaceColor',[0.8,0.8,0.8],...
      'facealpha',1);
hold on
% plot(dateVector,avgExp1*100,...
%        'color',[0.5 0.5 0.5],...
%        'linewidth',2)
plot(dateVector(B1),avgExp2(B1)*100,...
       'color',[0 0 0],...
       'linewidth',2)
plot(dateVector(B2),avgExp2(B2)*100,...
       'color',[0 0 0],...
       'linewidth',2,...
       'linestyle',':')
ylabel('Prepayments','interpreter','latex')

% Layout 
ax=gca;
ax.YAxis.TickLabelInterpreter='latex';
ax.YAxis.Label.Interpreter='latex';
ax.YAxis.TickLength=[0,0];
ax.YAxis.Label.FontSize=14;
ax.YAxis.FontSize=14;
ax.YAxis.Color=[0 0 0];
ax.XAxis.TickLabelInterpreter='latex';
ax.XAxis.Label.Interpreter='latex';
ax.XAxis.TickLength=[0,0];
ax.XAxis.Label.FontSize=14;
ax.XAxis.FontSize=14;

% Legend stuff
legend({'Prepayments\quad',...
        'Fitted Prepayment\quad',...
        'Out of sample'},...
        'interpreter','latex',...
        'location','northoutside',...
        'orientation','horizontal',...
        'fontsize',12,...
        'box','off')
    
% Fix dates for x ticks
xlim([dateVector(1)-90 dateVector(end)+90])
datetick;

%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\prepayOutOfSample')

% Residual plot
f=figure('color','w','position',[360   305   600   227]);
avgerror=nan(size(isins,1),1);
for i=1:size(isins,1)
    s=scatter(dateVector,(prepayments(:,i)-expPrepay1(:,i))*100);
    s.MarkerEdgeColor=ones(1,3)*0.5;
    s.Marker='.';
    hold on
end
residuals=100*reshape(prepayments-expPrepay1,size(isins,1)*size(dateVector,1),1);
allDates=repmat(dateVector,size(isins,1),1);
allDates(isnan(residuals))=[];
fprintf('%s%%',num2str(100*sum(allDates<=datenum('2011-01-01'))/size(allDates,1)))
residuals(isnan(residuals))=[];
residualsSort=sort(residuals);
N=size(residuals,1);
confBand=residualsSort(floor(N*[0.025 0.975]));
plot([dateVector(1),dateVector(end)],[confBand(1) confBand(1)],'color',ones(1,3)*0,'linestyle','--','linewidth',2)
plot([dateVector(1),dateVector(end)],[confBand(2) confBand(2)],'color',ones(1,3)*0,'linestyle','--','linewidth',2)
ylabel('Estimation Error')

% Draw average error
errors=(prepayments-expPrepay1)*100;
avgerror=nan(size(dateVector,1),1);
for i=1:size(dateVector,1)
    temp=errors(i,:);
    B=~isnan(temp);
    avgerror(i)=mean(temp(B));
end
%plot(dateVector,avgerror,'linewidth',2)

% Layout 
ax=gca;
ax.YAxis.TickLabelInterpreter='latex';
ax.YAxis.Label.Interpreter='latex';
ax.YAxis.TickLength=[0,0];
ax.YAxis.Label.FontSize=14;
ax.YAxis.FontSize=14;
ax.YAxis.Color=[0 0 0];
ax.XAxis.TickLabelInterpreter='latex';
ax.XAxis.Label.Interpreter='latex';
ax.XAxis.TickLength=[0,0];
ax.XAxis.Label.FontSize=14;
ax.XAxis.FontSize=14;
datetick
ax.Box='on';
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\prepayResiduals')

% Residual plot
f=figure('color','w','position',[360   305   600   227]);
for i=1:83
    s=scatter(dateVector,(marketProbs(:,i)-expPrepayProb1(:,i))*100);
    s.MarkerEdgeColor=ones(1,3)*0.5;
    s.Marker='.';
    
    hold on
end
residuals=100*reshape(marketProbs-expPrepayProb1,size(isins,1)*size(dateVector,1),1);
residuals(isnan(residuals))=[];
residualsSort=sort(residuals);
N=size(residuals,1);
confBand=residualsSort(floor(N*[0.025 0.975]));
plot([dateVector(1),dateVector(end)],[confBand(1) confBand(1)],'color',ones(1,3)*0,'linestyle','--','linewidth',2)
plot([dateVector(1),dateVector(end)],[confBand(2) confBand(2)],'color',ones(1,3)*0,'linestyle','--','linewidth',2)

% Layout 
ax=gca;
ax.YAxis.TickLabelInterpreter='latex';
ax.YAxis.Label.Interpreter='latex';
ax.YAxis.TickLength=[0,0];
ax.YAxis.Label.FontSize=14;
ax.YAxis.FontSize=14;
ax.YAxis.Color=[0 0 0];
ax.XAxis.TickLabelInterpreter='latex';
ax.XAxis.Label.Interpreter='latex';
ax.XAxis.TickLength=[0,0];
ax.XAxis.Label.FontSize=14;
ax.XAxis.FontSize=14;
datetick
ax.Box='on';

function val=GMMform2(theta,prepayments,XstarMatrix,LAMBDA,abovePar)
val=0;
B=~isnan(prepayments);
for i=1:size(XstarMatrix,3)
    val=val+GMMform(theta,prepayments(B(:,i),i),XstarMatrix(B(:,i),:,i),LAMBDA,abovePar(B(:,i),i));
end
end

function [val,expectedPrepayRates]=GMMform(theta,prepayments,Xstar,LAMBDA,abovePar)
% Get parameters
A=theta(1);
B=theta(2);
lambda=theta(3);
dt=0.25;

% Interpolate Xstar given lambda
XstarInt=nan(size(Xstar,1),1);
for i=1:size(Xstar,1)
    XstarInt(i)=hermiteInterpolation(LAMBDA,Xstar(i,:)',lambda);
end

% Calculate prepayments
expectedPrepayRates=modelPrepayments(A,B,lambda,XstarInt,1/dt);

% If below par then zero
expectedPrepayRates(~abovePar)=0;

% Calculate quadratic form
f=prepayments-expectedPrepayRates;
val=f'*f;
end