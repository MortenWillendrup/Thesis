%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Expected return table
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
scenarios=[-300 -200 -100 -50 -25 0 25 50 100 200 300]';
mat=nan(size(mkt.bonds.isins));
coupon=nan(size(mkt.bonds.isins));
prices=nan(size(mkt.bonds.isins));
outstanding=nan(size(mkt.bonds.isins));
OAS=nan(size(mkt.bonds.isins));
OAD=nan(size(mkt.bonds.isins));
OAC=nan(size(mkt.bonds.isins));
expectedReturn=nan(size(scenarios,1),size(mkt.bonds.isins,1));
accPrepay=nan(size(scenarios,1),size(mkt.bonds.isins,1));
for j=size(expectedReturn,2):-1:1
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

    % Outstanding larger than 500m and coupon less than 6
    if outstanding(j)>500000000&&coupon(j)<0.05

        % Get OAS
        OAS(j)=model.getOAS(shortRate,prices(j));
        [OAD(j),OAC(j)]=model.keyfigures(shortRate);%+OAS(j)/10000);

        % Expected return
        [expectedReturn(:,j),accPrepay(:,j)]=model.expReturn(shortRate,scenarios);%,OAS(j)/10000,0);
    end
end

% Generate bond names
names=cell(size(prices));
for i=1:size(prices,1)
    names{i}=sprintf('%s''%s',num2str(coupon(i)*100),right(datestr(mat(i),10),2));
end

% Pick out bonds
B1=outstanding>500000000&coupon<0.05;
uniqueNames=unique(names(B1));

% expectedReturn2
expectedReturn2=nan(size(expectedReturn,1),size(uniqueNames,1));
accPrepay2=nan(size(expectedReturn,1),size(uniqueNames,1));
for i=1:size(uniqueNames,1)
    B=strcmpi(names,uniqueNames(i))&~isnan(expectedReturn(1,:)');
    expectedReturn2(:,i)=mean(expectedReturn(:,B),2);
    accPrepay2(:,i)=mean(accPrepay(:,B),2);
end

% Plot
grid=(-300:300)';
selected=[1;2;4;11];
exp2plot=nan(size(grid,1),size(selected,1));
for i=1:size(selected,1)
    exp2plot(:,i)=hermiteInterpolation(scenarios,expectedReturn2(:,selected(i)),grid);
end

f=latexPlot('x',grid,...
    'y',exp2plot*100,...
    'legend',uniqueNames(selected),...
    'xlabel','Short rate shock (Bp.)',...
    'ylabel','Expected Return (\%)');
ylim([-15.5 15.5])
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\expectedReturns')