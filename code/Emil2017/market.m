classdef market < matlab.mixin.Copyable % Handle class automatically updates properties
    
    properties(Access=public)        
        bonds;
    end
    
    properties(Dependent)
        % Curves and Tenors
        ciborCurve=[];
        ciborTenors=[];
        swapCurve=[];
        swapTenors=[];
        citaCurve=[];
        citaTenors=[];
        capCurve=[];
        capTenors=[];
        marketTenors=[];
        
        % Date
        date=[];
        
        % Parameters
        kappa=[];
        sigma=[];
        
        % Shock
        shock={};
    end
    
    properties(Access=private)
        % Data structure
        data;
        parameters;
        dataFields;
        ciborB;
        swapB;
        citaB;
        capB;
        
        % Raw zero and forward rates
        zeroRates=[];
        zeroTenors=[];
        forwardRates=[];
        forwardDerivative=[];
        
        % Current Date
        currentDate=[];
        
        % Shock info matrix
        shockContainer={};
    end
    
    
    
    methods
        % Constructor
        function obj=market(varargin)
            
            % Input 
            for i=1:2:nargin
                switch varargin{i}
                    case 'date'
                        obj.date=datenum(varargin{i+1});
                end
            end
            
            % Load market data
            S=load('/users/helmig/dropbox/ku/speciale/Data/ratesContainer','data');
            obj.data=S.data;clear S;
            
            % Sort away illiquid swap rates and hardcode away any strage rates
            obj.data=rmfield(obj.data,'SWAP11Y');
            obj.data=rmfield(obj.data,'SWAP13Y');
            obj.data=rmfield(obj.data,'SWAP14Y');
            B=obj.data.SWAP9Y.dates==730919;
            obj.data.SWAP9Y.rates(B)=nan;
            
            % Fix 20Y and 30Y swaps
            d1=datenum('1998-07-27'); 
            d2=datenum('2002-03-01');
            B=obj.data.SWAP10Y.dates==d1;
            spread1=obj.data.SWAP20Y.rates(B)-obj.data.SWAP10Y.rates(B);
            spread2=obj.data.SWAP30Y.rates(B)-obj.data.SWAP10Y.rates(B);
            B=obj.data.SWAP10Y.dates==d2;
            spread3=obj.data.SWAP20Y.rates(B)-obj.data.SWAP10Y.rates(B);
            spread4=obj.data.SWAP30Y.rates(B)-obj.data.SWAP10Y.rates(B);
            B=(obj.data.SWAP30Y.dates>=d1)&(obj.data.SWAP30Y.dates<=d2);
            obj.data.SWAP20Y.rates(B)=obj.data.SWAP10Y.rates(B)+spread1+cumsum(ones(sum(B),1)*(spread3-spread1)/sum(B));
            obj.data.SWAP30Y.rates(B)=obj.data.SWAP10Y.rates(B)+spread2+cumsum(ones(sum(B),1)*(spread4-spread2)/sum(B));
            
            % Sort away redundant cibor rates
            obj.data=rmfield(obj.data,'CIBOR07M');
            obj.data=rmfield(obj.data,'CIBOR08M');
            obj.data=rmfield(obj.data,'CIBOR09M');
            obj.data=rmfield(obj.data,'CIBOR10M');
            obj.data=rmfield(obj.data,'CIBOR11M');
            obj.data=rmfield(obj.data,'CIBOR12M');
            
            % Save fields
            obj.dataFields=fields(obj.data);
            
            % Create boolean to find cibor, swap, cita and caps
            obj.ciborB=~cellfun(@isempty,(strfind(lower(obj.dataFields),'cibor')));
            obj.swapB=~cellfun(@isempty,(strfind(lower(obj.dataFields),'swap')));
            obj.citaB=~cellfun(@isempty,(strfind(lower(obj.dataFields),'cita')));
            obj.capB=~cellfun(@isempty,(strfind(lower(obj.dataFields),'cap')));
            
            % Load calibrated and estimated parameters 
            S=load('/users/helmig/dropbox/ku/speciale/Data/parameters','data');
            obj.parameters=S.data;
            
%             % Load bonds
%             S=load('/users/helmig/dropbox/ku/speciale/Data/bonds','bonds');
%             obj.bonds=S.bonds;
            
            S=load('/users/helmig/dropbox/ku/speciale/Data/dataFromSEB');
            bonds=struct();
            bonds.isins=S.dataFromSEB.bonds(:,1);
            vpids=str2double(S.dataFromSEB.bonds(:,6));
            prepayDates=S.dataFromSEB.datesPrepay;
            prepayRates=S.dataFromSEB.prepayMatrix;
            isins2remove=false(size(bonds.isins));
            for i=1:size(bonds.isins,1)
                % Prices
                name=bonds.isins{i};
                tempPrices=S.dataFromSEB.priceMatrix(:,i);
                B=~isnan(tempPrices);
                bonds.(name).prices=tempPrices(B);
                bonds.(name).datesPrices=S.dataFromSEB.datesPrices(B);
                
                % Outstanding
                tempOA=fillmissing(S.dataFromSEB.outstanding(:,i),'previous');
                bonds.(name).outstanding=tempOA(B);               
                
                % Prepayments
                B=vpids==str2double(S.dataFromSEB.bonds(i,6));
                temp=prepayRates(:,B);
                B=~isnan(temp);
                bonds.(name).prepayments=temp(B);
                temp=datevec(prepayDates(B));                
                bonds.(name).datesPrepayments=prepayDates(B)-temp(:,3)+1;
                
                % Bond info
                bonds.(name).bank=S.dataFromSEB.bonds{i,8};
                bonds.(name).maturity=S.dataFromSEB.bonds{i,4};                
                bonds.(name).coupon=S.dataFromSEB.bonds{i,5};  
                
                % Select empty ones
                if isempty(bonds.(name).prices)
                    isins2remove(i)=true;
                    bonds=rmfield(bonds,name);
                end                
            end
            % Remove empty ones
            bonds.isins(isins2remove)=[];   
            
            % Also remove DK0009753113
            bonds.('DK0009753113')=[];
            bonds.isins(strcmpi(bonds.isins,'DK0009753113'))=[];
            
            % Save
            obj.bonds=bonds;
        end
        
        function [dateVector,values]=timeSeries(obj,name,startDate,endDate)
            % Find relevant rate
            B=strcmpi(obj.dataFields,name);
            
            % ... or price
            B2=strcmpi(obj.bonds.isins,name);
            
            % Datenumer
            startDate=datenum(startDate);
            endDate=datenum(endDate);
            
            % Create datevector
            dateVector=(startDate:endDate)';
            values=nan(size(dateVector));
            
            % Create data series
            if any(B)
                % Get data
                dates=obj.data.(obj.dataFields{B}).dates;
                [B2,IX]=ismember(dateVector,dates);
                values(B2)=obj.data.(obj.dataFields{B}).rates(IX(B2))./100;
            elseif any(B2)
                dates=obj.bonds.(name).datesPrices;
                [B2,IX]=ismember(dateVector,dates);
                values(B2)=obj.bonds.(name).prices(IX(B2));
            else
                dateVector=[];
                values=[];
            end
            
        end
        
        function curve=zeroCurve(obj,tenors)
            % If curve has not been calibrated then calibrate
            if isempty(obj.zeroRates)
                obj.swapCalibration;
            end
            
            % With the calibrated curve in place we provide the 
            % interpolated curve
            curve=hermiteInterpolationFwd(obj.zeroTenors,...
                                          obj.zeroRates,...
                                          tenors);
        end
        
        function curve=forwardCurve(obj,tenors)
            % If curve has not been calibrated then calibrate
            if isempty(obj.zeroRates)
                obj.swapCalibration;
            end
            
            % With the calibrated curve in place we provide the 
            % interpolated curve
            [~,curve]=hermiteInterpolationFwd(obj.zeroTenors,...
                                                obj.zeroRates,...
                                                tenors);
        end
        
        function curve=forwardCurveDiff(obj,tenors)
            % If curve has not been calibrated then calibrate
            if isempty(obj.zeroRates)
                obj.swapCalibration;
            end
            
            % With the calibrated curve in place we provide the 
            % interpolated curve. We interpolate the forward curve in order
            % to get a more smooth derivative
            [~,~,~,curve]=hermiteInterpolation(obj.zeroTenors,...
                                                obj.forwardRates,...
                                                tenors);
        end
        
        
    end
    
    % Methods for dependent properties. 
    methods 
        % CIBOR Curve
        function curve=get.ciborCurve(obj)
            curve=getCurve(obj,'cibor')./100;
        end
        
        % CIBOR Tenors
        function tenors=get.ciborTenors(obj)
            tenors=getTenors(obj,'cibor');
        end 
         
        % SWAP Curve
        function curve=get.swapCurve(obj)
            curve=getCurve(obj,'swap')./100;
        end
        
        % SWAP Tenors
        function tenors=get.swapTenors(obj)
            tenors=getTenors(obj,'swap');
        end 
         
        % CITA Curve
        function curve=get.citaCurve(obj)
            curve=getCurve(obj,'cita')./100;
        end
        
        % CITA Tenors
        function tenors=get.citaTenors(obj)
            tenors=getTenors(obj,'cita');
        end 
         
        % CAP Curve
        function curve=get.capCurve(obj)
            curve=getCurve(obj,'cap');
        end
        
        % CAP Tenors
        function tenors=get.capTenors(obj)
            tenors=getTenors(obj,'cap');
        end
        
        % Market Tenors
        function tenors=get.marketTenors(obj)
            tenors=[obj.ciborTenors;obj.swapTenors];
        end
        
        % Date get and set functions
        function val=get.date(obj)
            val=obj.currentDate;
        end
        function set.date(obj,val)
            % Set date
            obj.currentDate=datenum(val);
            
            % Clear any calibrated curves
            obj.zeroRates=[];
            obj.forwardRates=[];
        end
        
        % Shocks
        function val=get.shock(obj)
            val=obj.shockContainer;
        end
        function set.shock(obj,val)
            % Set date
            obj.shockContainer=val;
            
            % Clear any calibrated curves
            obj.zeroRates=[];
            obj.forwardRates=[];
        end
        
        % Get kappa
        function val=get.kappa(obj)
            if obj.date<obj.parameters.dateVector(1)
                val=mean(obj.parameters.kappa);return; %Extrapolate 
            else
                B=obj.parameters.dateVector==obj.date;
            end
            if any(B)
                val=obj.parameters.kappa(B);
            else
                val=nan;
                warning('No parameters available for date %s',datestr(obj.date))
            end
        end
        
        % Get sigma
        function val=get.sigma(obj)
            if obj.date<obj.parameters.dateVector(1)
                val=mean(obj.parameters.sigma);return; %Extrapolate 
            else
                B=obj.parameters.dateVector==obj.date;
            end
            if any(B)
                val=obj.parameters.sigma(B);
            else
                val=nan;
                warning('No parameters available for date %s',datestr(obj.date))
            end
        end
        
    end
    
    methods(Access=private)
        function curve=getCurve(obj,curveName)
            switch lower(curveName)
                case 'swap'
                    B=obj.swapB;
                case 'cibor'
                    B=obj.ciborB;
                case 'cita'
                    B=obj.citaB;
                case 'cap'
                    B=obj.capB;
            end
            
            % Get fields
            fields=obj.dataFields(B);
            
            % Preallocate 
            curve=nan(size(fields));
            tenors=nan(size(fields));
            
            % If curve is not available (holiday/weekend) then take day
            % before
            k=0;
            while sum(~isnan(curve))<5
                % Loop to get data
                for i=1:size(fields,1)
                    B=obj.data.(fields{i}).dates==(obj.currentDate-k);
                    curve(i)=obj.data.(fields{i}).rates(B);
                    tenors(i)=obj.data.(fields{i}).tenor;
                end
                k=k+1;
            end
            
            % Sort by tenor
            [~,IX]=sort(tenors);
            curve=curve(IX); 
            
            % Remove any nan
            B=isnan(curve);
            curve(B)=[];
            
            % Add any shock
            if ~isempty(obj.shock)
                fields=fields(IX);fields(B)=[];
                for i=1:size(obj.shock,1)
                    B=strcmpi(fields,obj.shock{i,1});
                    curve(B)=curve(B)+obj.shock{i,2}*100;
                end
            end
            
        end
    
        function tenors=getTenors(obj,curveName)
            switch lower(curveName)
                case 'swap'
                    B=obj.swapB;
                case 'cibor'
                    B=obj.ciborB;
                case 'cita'
                    B=obj.citaB;
                case 'cap'
                    B=obj.capB;
            end
            
            % Get fields
            fields=obj.dataFields(B);
            
            % Preallocate 
            curve=nan(size(fields));
            tenors=nan(size(fields));
            
            % If curve is not available (holiday/weekend) then take day
            % before
            k=0;
            while sum(~isnan(curve))<5
                % Loop to get data
                for i=1:size(fields,1)
                    B=obj.data.(fields{i}).dates==(obj.currentDate-k);
                    curve(i)=obj.data.(fields{i}).rates(B);
                    tenors(i)=obj.data.(fields{i}).tenor;
                end
                k=k+1;
            end
            
            % Sort by tenor
            [tenors,IX]=sort(tenors);
            curve=curve(IX); 
            
            % Remove any nan
            B=isnan(curve);
            tenors(B)=[]; 
        end
        
        % Function for calibrating swap curve
        function obj=swapCalibration(obj)

            % Create zero tenors
            obj.zeroTenors=[obj.ciborTenors;obj.swapTenors];

            % CIBOR is straight forward y=-1/tau*log(1/(1+delta*L))
            ciborZero=-1./obj.ciborTenors.*log(1./(1+obj.ciborTenors.*obj.ciborCurve));

            % Get number of swap tenors
            N=size(obj.swapTenors,1);

            % Initial guess is just set to the swap rates
            guess=obj.swapCurve;
            swapFit=obj.swapCurve+1;

            % Find all relevant tenors
            tenors=(0.5:0.5:max(obj.swapTenors))';
            while sum(abs(obj.swapCurve-swapFit))*10000>1 % Max 0.1 BPS deviation
                % Interpolate zero rates
                zero=hermiteInterpolationFwd(obj.zeroTenors,[ciborZero;guess],tenors);

                % Calculate ZCBs
                ZCB=exp(-tenors.*zero);

                % Calculate new zero rates through Hagan & West algorithm 
                sumP=nan(N,1);
                for i=1:N
                    B=tenors<obj.swapTenors(i);
                    sumP(i)=sum(ZCB(B));
                end
                temp=(1-obj.swapCurve.*sumP*0.5)./(1+obj.swapCurve*0.5);

                % Compute new guess
                guess=-1./obj.swapTenors.*log(temp);

                % Calculate theoretical swap rates
                swapFit=nan(N,1);
                for i=1:N
                    T=obj.swapTenors(i);
                    B=tenors==T;
                    numerator=1-ZCB(B);
                    B=ismember(tenors,0.5:0.5:T);
                    denominator=0.5*sum(ZCB(B));
                    swapFit(i)=numerator/denominator;
                end
            end

            % Compute output
            [obj.zeroRates,obj.forwardRates]=...
                hermiteInterpolationFwd(obj.zeroTenors,...
                                       [ciborZero;guess],...
                                       obj.zeroTenors);

        end
        
    end
    
end

