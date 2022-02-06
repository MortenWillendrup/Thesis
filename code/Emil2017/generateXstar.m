%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Generate X-star over time for multiple bonds
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
mkt=market('date','2017-06-01');

% Load isins
isins=mkt.bonds.isins;
N3=size(isins,1);

% Load dates
minDate=mkt.bonds.(isins{1}).datesPrepayments(1);
for i=2:N3
    temp=mkt.bonds.(isins{i}).datesPrepayments;
    if ~isempty(temp)
        minDate=min(minDate,temp(1));
    end
end
maxDate=mkt.bonds.(isins{1}).datesPrepayments(end);
for i=2:N3
    temp=mkt.bonds.(isins{i}).datesPrepayments;
    if ~isempty(temp)
        maxDate=max(maxDate,temp(end));
    end
end
dateVector=(minDate:maxDate)';
dateInfo=datevec(dateVector);
B=mod(dateInfo(:,2)-1,3)==0&dateInfo(:,3)==1;
dateVector(~B)=[];
N=size(dateVector,1);

% Load maturity and coupon info
info=nan(N3,2);
for i=1:N3
    info(i,1:2)=[datenum(mkt.bonds.(isins{i}).maturity),mkt.bonds.(isins{i}).coupon];
end

% Some info are duplicates. Indicate the first one in order to not
% calculate stuff twice or more
B=false(N3,1);
for i=1:N3
    if sum((info(i,1)==info(1:i,1))&(info(i,2)==info(1:i,2)))==1
        B(i)=true;
    end
end
info(~B,:)=[];
N3=size(info,1);

% Load 1M Cibor
[x,y]=mkt.timeSeries('cibor01m',dateVector(1)-10,dateVector(end)+10);
y=fillmissing(y,'previous');
cibor=nan(size(dateVector));
for i=1:size(dateVector,1)
    B=x==dateVector(i);
    cibor(i)=y(B);
end

% Parameters
X=0.15;
F0=100;
n=4;
T=30;
lambda=0.5202;

% Initiate model
model=cStantonHW('X',X,...
                 'F0',F0,...
                 'n',n,...
                 'T',T,...
                 'lambda',lambda,...
                 'market',mkt);


% Calculate Xstar over time for different combinations of coupons and lambda
LAMBDA=(0:0.1:4)';
N2=size(LAMBDA,1);
Xstar=nan(N,N2,N3);
for i=1:N
    % Setting up model
    fprintf('Setting up model for date %s of %s\n',num2str(i),num2str(N))
    
    % Update market 
    mkt.date=dateVector(i);
    model.market=mkt;
    
    % Loop over bonds 
    for k=1:N3
        
        % Calculate 
        fprintf('Calculating X for bond %s of %s\n',num2str(k),num2str(N3))
        
        % Set model parameters
        model.maturity=round((info(k,1)-dateVector(i))/365*4,0)/4;
        model.R=info(k,2);
        
        % Prepare model for pricing
        model.prepare;

        % Loop over lambda
        for j=1:N2
            % Set lambda
            model.lambda=LAMBDA(j);

            % Calculate Xstar given lambda
            lb=0;
            ub=1;    
            l=999;
            counter=0;
            model.X=0;
            while abs(l-100*(1+model.X))>0.01&&counter<14
                mb=0.5*(lb+ub);    
                model.X=mb;
                [Ma,Ml,r]=model.pricingFD;
                l=hermiteInterpolation(r,Ml,cibor(i));
                if l>100*(1+model.X)
                    lb=mb;
                else
                    ub=mb;
                end
                counter=counter+1;

                % Terminate while loop if we hit 0 og 1
                if round(mb,3)==0
                    mb=0;break;
                elseif round(mb,3)==1
                    mb=1;break;
                end
            end

            % Save result
            Xstar(i,j,k)=mb;
        end
    end
    % Clear screen
    clc;
end

save('/users/helmig/dropbox/ku/speciale/Data/estimationData','dateVector','Xstar','info','LAMBDA')
