%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Expected return figures
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
monitorDate=datenum('2017-01-01');
mkt=market('date',monitorDate);

% Load short rate
[x,y]=mkt.timeSeries('cibor06M',monitorDate-10,monitorDate+10);
y=fillmissing(y,'previous');
shortRate=y(x==monitorDate);

% Parameters
F0=100;
n=4;
lambda=0.5773;
Alpha=0.0612;
Beta=2.2511;

% Initiate model
model=cStantonHW('F0',F0,...
                 'n',n,...
                 'lambda',lambda,...
                 'A',Alpha,...
                 'B',Beta,...
                 'market',mkt,...
                 'timesteps',3,...
                 'spacesteps',400,...
                 'rannacher',true,...
                 'smoothing',true);

% Load estimation data
data=load('/users/helmig/dropbox/ku/speciale/Data/estimationData');
info=data.info;
datesXstar=data.dateVector;
XstarCube=data.Xstar;
LAMBDA=data.LAMBDA;

% Loop over bonds
expectedReturn=nan(size(mkt.bonds.isins));
OAD=expectedReturn;
OAC=expectedReturn;
OAS=expectedReturn;
mat=expectedReturn;
coupon=expectedReturn;
prices=expectedReturn;
outstanding=expectedReturn;
rGrid=(-0.1:0.001:0.2)';
priceMatrix=nan(size(mkt.bonds.isins,1),size(rGrid,1));
for j=1:size(expectedReturn,1)
    % Load prepaymetns
    isin=mkt.bonds.isins{j};
    dates=mkt.bonds.(isin).datesPrepayments(1:end-3);
    mat(j)=datenum(mkt.bonds.(isin).maturity);

    % Remove dates in issuing period
    B=dates<=floor(mat(j)-30*365.25);
    dates(B)=[];
    
    % Pick out bond
    B=info(:,2)==mkt.bonds.(isin).coupon&info(:,1)==datenum(mkt.bonds.(isin).maturity);
    if any(B)
        XstarMat=XstarCube(:,:,B);
        Xstar=nan(size(XstarMat,1),1);
        for i=1:size(Xstar,1)
            Xstar(i)=hermiteInterpolation(LAMBDA,XstarMat(i,:)',lambda);
        end

        % Crop Xstar to prepayment dates
        B=ismember(datesXstar,dates);
        Xstar(~B)=[];	
        model.Xstar=Xstar;
    
    else
        model.Xstar=[];    
    end
    
    % Set coupon
    coupon(j)=mkt.bonds.(isin).coupon;
    model.R=coupon(j);
    %if ~(coupon(j)==0.015);continue;end;
    % Set maturity
    model.maturity=round(4*(mat(j)-monitorDate)/365.25,0)/4;

    % Get market price
    temp1=mkt.bonds.(isin).prices;
    temp2=mkt.bonds.(isin).outstanding;
    datesPrices=mkt.bonds.(isin).datesPrices;
    IX=sum(datesPrices<=monitorDate);
    if IX==0;continue;end
    prices(j)=temp1(IX);
    outstanding(j)=temp2(IX);

    % Get OAS
    OAS(j)=model.getOAS(shortRate,prices(j));
    [OAD(j),OAC(j)]=model.keyfigures(shortRate);%+OAS(j)/10000);

    % Expected return
    expectedReturn(j)=model.expReturn(shortRate,0);%,OAS(j)/10000,0);
    
    % Get price function
    %priceMatrix(j,:)=model.getPrice(rGrid)';
end

% Outstanding larger than 500m and coupon less than 6
B1=outstanding>5e8&coupon<=0.06;

% Calc VaR and ES
names=cell(size(OAD));
for i=1:size(OAD,1)
    names{i}=sprintf('%s''%s',num2str(coupon(i)*100),right(datestr(mat(i),10),2));
end

% Average over names
uniqueNames=unique(names(B1));
avgReturn=nan(size(uniqueNames));
avgOAD=nan(size(uniqueNames));
avgOAS=nan(size(uniqueNames));
for i=1:size(uniqueNames,1)
    B=strcmpi(names,uniqueNames(i));
    avgReturn(i)=mean(expectedReturn(B&B1));
    avgOAD(i)=mean(OAD(B&B1));
    avgOAS(i)=mean(OAS(B&B1));
end

% Ols
XX=[ones(size(avgOAS)),avgOAS*10000];
YY=100*avgReturn;
BETA=(XX'*XX)\(XX'*YY);
residuals=YY-XX*BETA;

% Fig
f=figure('color','w','position',[360   305   600   350]);
p=plot((40:360)',BETA(1)+BETA(2)*(40:360)');
p.Color=[0 0 0];
p.LineStyle='--';
hold on
scatter(avgOAS*10000,avgReturn*100,'MarkerEdgeColor','none','MarkerFaceColor',[0 0 0])
for i=1:size(avgOAS,1)    
    t=text(avgOAS(i)*10000+5,avgReturn(i)*100+0*0.1,...
        uniqueNames{i},...
        'fontsize',12,...
        'interpreter','latex',...
        'HorizontalAlignment','left');
end
%ylim([1.6 3])
%xlim([6 14])
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
xlabel('OAS (Bp.)')
ylabel('Expected Return (\%)')
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\OASfig3')



f=figure;
[~,ix]=sort(OAC,'descend');
randomP2=randomP(ix,:);
for i=1:size(OAD,1)
    histogram(randomP(i,:))
    xlim([71 137])
    pause(0.5)
    cla(f.Children)
end

% Figure
model.Xstar=[];
model.maturity=30;
model.R=0.015;
oas=model.getOAS(shortRate,96.82);
model.expReturn(shortRate,-200);

% Simulate return distributions
vol=sqrt(model.sigma^2/(2*model.kappa)*(1-exp(-2*model.kappa)))*100;

X=randn(10000,1)*vol;
B=X<0;
Y1=-2*X-1/2*0*X.^2;
Y2=-0.5*X-1/2*2.7*X.^2;
Y3=-5*X-1/2*5*X.^2;
subplot(1,3,1)
delta=(ceil(max(Y1))-floor(min(Y1)))/50;
binrng=floor(min(Y1)):delta:ceil(max(Y1));
counts1 = histc(Y1(B), binrng); 
counts2 = histc(Y1(~B), binrng); 
counts3 = counts1 + counts2; 
figure(1)
bar(binrng, counts3, 'b')
hold on
bar(binrng, counts1, 'y')
hold off
legend('Negative', 'Positive')
subplot(1,3,2)
delta=(ceil(max(Y2))-floor(min(Y2)))/50;
binrng=floor(min(Y2)):delta:ceil(max(Y2));
counts1 = histc(Y2(B), binrng); 
counts2 = histc(Y2(~B), binrng); 
counts3 = counts1 + counts2; 
figure(1)
bar(binrng, counts3, 'b')
hold on
bar(binrng, counts1, 'y')
hold off
legend('Negative', 'Positive')
subplot(1,3,3)
delta=(ceil(max(Y3))-floor(min(Y3)))/50;
binrng=floor(min(Y3)):delta:ceil(max(Y3));
counts1 = histc(Y3(B), binrng); 
counts2 = histc(Y3(~B), binrng); 
counts3 = counts1 + counts2; 
figure(1)
bar(binrng, counts3, 'b')
hold on
bar(binrng, counts1, 'y')
hold off
legend('Negative', 'Positive')



