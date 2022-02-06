%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Correlation study
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Get market
mkt=market;

% Interest rates and prices
startDate=datenum('1994-01-01');
endDate=datenum('2017-06-01');
rateNames={'Cibor01M'
    'Cibor06M'
    'Swap1Y'
    'Swap5Y'
    'Swap10Y'
    'DK0009334575'};

dateVector=(startDate:endDate)';
N1=size(dateVector,1);
N2=size(rateNames,1);
data=nan(N1,N2);
for i=1:N2
    [~,data(:,i)]=mkt.timeSeries(rateNames{i},dateVector(1),dateVector(end));
end

% Get zero rates
dateVector;
terms=[1;2;3;5;10;15;20;30];
yields=nan(size(dateVector,1),size(terms,1));
forwards=nan(size(dateVector,1),size(terms,1));
h=waitbar(0);
for i=1:size(dateVector,1)
    mkt.date=dateVector(i);
    yields(i,:)=mkt.zeroCurve(terms)';
    forwards(i,:)=mkt.forwardCurve(terms)';
    waitbar(i/size(dateVector,1),h)
end

% Weekly Diff
data=fillmissing(data,'previous');
B=any(isnan(data),2);
data(B,:)=[];
diffRates=data(15:14:end,:)-data(1:14:end-14,:);
sig=cov(diffRates(:,1:5));


cumsum(sort(eig(cov(diff(yields(1:30:end,:))))./sum(eig(cov(diff(yields(1:30:end,:))))),'descend'))
[eigVal,~,eigVec]=eig(cov(diff(yields(1:30:end,:)*100)));
