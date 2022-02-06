%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Lehman graphs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Get market
startDate=datenum('2007-01-01');
endDate=datenum('2010-01-01');
mkt=market('date',startDate);

% Get top N largest bonds around Lehman
N=3;
isins=mkt.bonds.isins;
outstanding=nan(size(isins));
for i=1:size(isins,1)
    dates=mkt.bonds.(isins{i}).datesPrices;
    B=dates<=startDate;
    if any(B)
        outstanding(i)=mkt.bonds.(isins{i}).outstanding(sum(B));
    end
end
B=isnan(outstanding);
isins(B)=[];
outstanding(B)=[];
[outstanding,ix]=sort(outstanding,'descend');
isins=isins(ix);
outstanding=outstanding(1:N);
isins=isins(1:N);

% Get prices
dateVector=(startDate:endDate)';
prices=nan(size(dateVector,1),N);
for i=1:N
    tempDates=mkt.bonds.(isins{i}).datesPrices;
    tempPrices=mkt.bonds.(isins{i}).prices;
    for j=1:size(dateVector,1)
        B=tempDates==dateVector(j);
        if any(B)
            prices(j,i)=tempPrices(B);
        end
    end
    prices(:,i)=fillmissing(prices(:,i),'previous');
end


