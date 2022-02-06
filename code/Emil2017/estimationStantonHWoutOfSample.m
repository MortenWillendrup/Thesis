%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Stanton estimation under Hull & White dynamics
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
mkt=market('date','2017-06-01');

% Load prepaymetns
isin=mkt.bonds.isins{1};
prepay=mkt.bonds.(isin).Extraordinary(1:end-1);
outstanding=mkt.bonds.(isin).Outstanding(1:end-1);
dates=mkt.bonds.(isin).Dates(1:end-1);
mat=mkt.bonds.(isin).Maturity;
time=round((mat-dates)/365*4,0)/4;
N=size(time,1);

% Load 1M Cibor
[x,y]=mkt.timeSeries('cibor01m',dates(1)-10,dates(end));
y=fillmissing(y,'previous');
cibor=nan(size(dates));
for i=1:N
    B=x==dates(i);
    cibor(i)=y(B);
end

% Load 10Y Swap
[x,y]=mkt.timeSeries('swap10Y',dates(1)-10,dates(end));
y=fillmissing(y,'previous');
swap10y=nan(size(dates));
for i=1:N
    B=x==dates(i);
    swap10y(i)=y(B);
end

% Parameters
X=0.15;
F0=100;
R=mkt.bonds.(isin).Coupon;
n=4;
T=30;
lambda=0.5202;

% Initiate model
model=cStantonHW('X',X,...
                 'F0',F0,...
                 'R',R,...
                 'n',n,...
                 'T',T,...
                 'lambda',lambda,...
                 'market',mkt);

% Generate model for each point in time
modelContainer=cell(N,1);
for i=1:N
    % Setting up model
    fprintf('Setting up model for date %s of %s',num2str(i),num2str(N))
    
    % Update market 
    mkt.date=dates(i);
    model.market=copy(mkt);
    
    % Update time to maturity
    model.T=time(i);
    
    % Prepare model for pricing
    model.prepare;
    
    % Save model    
    modelContainer{i}=copy(model);
    clc;
end


% Calculate X^* over time and over lambda
LAMBDA=(0:0.1:4)';
N2=size(LAMBDA,1);
Xstar=nan(N,N2);
for i=1:N
    % Setting up model
    fprintf('Calculating X for date %s of %s\n',num2str(i),num2str(N))
    model=modelContainer{i};
    
    % Loop over lambda
    for j=1:N2
        % Set lambda
        model.lambda=LAMBDA(j);
        
        % Calculate Xstar given lambda
        if j==1
            lb=0;
            ub=1;    
        else
            % If we allready calculated calculated Xstar once then we know
            % that Xstar is decreasing in lambda so we set the upper bound
            % to the previous Xstar
            lb=0;
            ub=mb;
        end
        l=101;
        k=0;
        while abs(l-100*(1+model.X))>0.01&&k<14
            mb=0.5*(lb+ub);    
            model.X=mb;
            [Ma,Ml,r]=model.pricingFD;
            l=hermiteInterpolation(r,Ml,cibor(i));
            if l>100*(1+model.X)
                lb=mb;
            else
                ub=mb;
            end
            k=k+1;
            
            % Terminate while loop if we hit 0 og 1
            if round(mb,3)==0
                mb=0;break;
            elseif round(mb,3)==1
                mb=1;break;
            end
        end
        
        % Save result
        Xstar(i,j)=mb;
    end
    % Clear screen
    clc;
end

[swapDates,swap]=mkt.timeSeries('swap10y',dates(1),dates(end));
% 
% yyaxis left
% plot(dates,Xstar);datetick;
% yyaxis right
% plot(dates,prepay)

% Estimate
THETA=nan(size(Xstar,1),3);
expectedPrepay=nan(size(Xstar,1),1);
for i=29:size(Xstar,1)
    A=2;
    B=10;
    theta=[A;B;2];
    lb=zeros(3,1);
    ub=[100;100;3];
    tic
    theta=fmincon(@(theta)GMMform(theta,prepay(1:i),Xstar(1:i,:),LAMBDA),theta,[],[],[],[],lb,ub)
    toc
    [val,temp]=GMMform(theta,prepay(1:i),Xstar(1:i,:),LAMBDA);
    THETA(i,:)=theta';
    expectedPrepay(i)=temp(end);
end
% Plot prepayments
f=figure('color','w','position',[360   305   600   273]);

% Plot prepayments 
yyaxis left;
b=bar(dates,prepay*100,...
      'barWidth',1,...
      'FaceColor',[0.8,0.8,0.8],...
      'facealpha',1);
hold on
plot(dates,expectedPrepay*100,...
       'color',[0 0 0],...
       'linewidth',2)
ylabel('Prepayments','interpreter','latex')

% Plot swap
yyaxis right;
p=plot(swapDates(1:7:end),swap(1:7:end)*100,...
       'color',[0 0 0]);
ylabel('10Y Swap (\%)','interpreter','latex')
xlabel('Year','interpreter','latex')

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
        'Estimated Prepayment (Left)\quad',...
        '10Y Swap (Right)'},...
        'interpreter','latex',...
        'location','northoutside',...
        'orientation','horizontal',...
        'fontsize',12,...
        'box','off')
    
% Fix dates for x ticks
dateVector=(dates(1):dates(end))';
dateMatrix=datevec(dateVector);
B=dateMatrix(:,2)==1&dateMatrix(:,3)==1;
xDateTicks=dateVector(B);
xDateTicks=xDateTicks(1:2:end);
xDateString=num2str(datestr(xDateTicks,10));
ax.XAxis.TickValues=xDateTicks;
ax.XAxis.TickLabels=xDateString;
xlim([dates(1)-90 dates(end)+90])

%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\prepayDK0009753469')

%%% Borrower distribution over time

% Get parameters
A=theta(1);
B=theta(2);
lambda=theta(3);

% Interpolate Xstar given lambda
XstarInt=nan(size(Xstar,1),1);
for i=1:size(Xstar,1)
    XstarInt(i)=hermiteInterpolation(LAMBDA,Xstar(i,:)',lambda);
end

% Selected dates for which we want distribution
distDates={'2004-01-01';'2008-01-01';'2012-01-01';'2016-01-01'};
xGrid=(0:0.01:1)';
PDF=nan(size(xGrid,1),size(distDates,1));
CDF=nan(size(xGrid,1),size(distDates,1));
for j=1:size(distDates,1)
    b=dates<=datenum(distDates(j));
    for i=1:size(xGrid,1)
        [PDF(i,j),CDF(i,j)]=borrowerDistribution(xGrid(i),XstarInt(b),lambda,A,B,n);
    end
end
f=latexPlot('x',xGrid*100,'y',PDF,...
            'xlabel','Prepayment Costs, $X$ (\%)',...
            'ylabel','Density function, $f_t^X$',...
            'legend',strcat('\textit{',strleft(distDates,4),'}$\quad$'),...
            'location','northoutside',...
            'orientation','horizontal');
ylim([0 5])
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\prepayDistDK0009753469')

%%% Calculate bond prices
modelPrices=nan(N,1);
for i=1:N
    % Setting up model
    fprintf('Calculating prices for date %s of %s\n',num2str(i),num2str(N))
    model=modelContainer{i};
    
    % Set model parameters
    model.lambda=lambda;
    
    % Calculate borrower distribution 
    xGrid=(0:0.1:1)';
    CDF=nan(size(xGrid));
    for j=1:size(xGrid,1)
        [~,CDF(j)]=borrowerDistribution(xGrid(j),XstarInt(1:i),lambda,A,B,n);
    end
    
    % Compute probability vector
    probs=CDF(2:end)-CDF(1:end-1);
    xGrid=(xGrid(2:end)+xGrid(1:end-1))./2;
    
    % Get prices
    prices=nan(size(xGrid));
    for j=1:size(xGrid,1)
        model.X=xGrid(j);
        [Ma,~,r]=model.pricingFD;
        prices(j)=hermiteInterpolation(r,Ma,cibor(i));
    end
    
    % Weightes prices
    modelPrices(i)=prices'*probs;
    
    % Clear screen
    clc;
end

function [val,expectedPrepayRates]=GMMform(theta,prepayments,Xstar,LAMBDA)
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

% Calculate prepayment probabilities
prepayProbs=nan(size(prepayments));
for i=1:size(XstarInt,1)
    prepayProbs(i)=prepaymentProbability(XstarInt(i),XstarInt(1:i-1),lambda,A,B,1/dt);
end

% Calculate prepayment rates
remaining=1-cumsum(prepayProbs);
expectedPrepayRates=prepayProbs./remaining;

% Calculate quadratic form
f=prepayments-expectedPrepayRates;
val=f'*f;
end