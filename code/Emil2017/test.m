%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Structural changes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
mkt=market;
isin='DK0009753469';
mat=datenum(mkt.bonds.(isin).maturity);
R=mkt.bonds.(isin).coupon;

% Load bond prices
dateVector=mkt.bonds.(isin).datesPrices;
bondPrices=mkt.bonds.(isin).prices;

% Load prepayments
datePre=mkt.bonds.(isin).datesPrepayments;
prepayments=mkt.bonds.(isin).prepayments;

% Load estimation data
S=load('/users/helmig/dropbox/ku/speciale/Data/estimationData','dateVector','Xstar','info','LAMBDA');

% Pick out Xstar
B=S.info(:,1)==mat&S.info(:,2)==R;
XstarMat=S.Xstar(:,:,B);
dates=S.dateVector;

% Hermite
Xstar=nan(size(XstarMat,1),1);
for i=1:size(Xstar,1)
    Xstar(i)=hermiteInterpolation(S.LAMBDA,XstarMat(i,:)',0.5);
end

% Correct Xstar
for i=1:size(Xstar,1)
    B=dateVector<=dates(i);
    if any(B)
        price=bondPrices(sum(B));
        if price<100
            Xstar(i)=0;
        end
    end
end

% Plot
bar(datePre,prepayments);hold on
plot(dates,modelPrepayments(0.2,4,0.5,Xstar,4)*100)
