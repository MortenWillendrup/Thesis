classdef cStantonHW < matlab.mixin.Copyable
    % Pricing in the Stanton (1995) model
    
    properties (Access = public)
        % Preallocate non-dependent model parameters
        X=[];
        F0=[];
        n=[];
        T=[];
        lambda=[];
        Alpha=[];
        Beta=[];
        delivery=true;
        feeRateShort=0.00;
        feeRateLong=0.00;
    end
    
    properties (Access = private)
        % FD matrices
        matrixContainer={};
        
        % Private Model parameters
        thetaPrivate=[];
        XstarPrivate=[];
        RPrivate=[];
        OASPrivate=[];
        
        % Private market
        mkt=[];
        
        % Preallocate finite difference settings
        % Set discretisation parameters
        n2=1; % Time steps between payment dates
        t0=0;
        r_min=-0.1;
        r_max=0.3;
        N=200;
        RannacherTimeStepping=true;
        smoothPayoff=true;
        RannacherSteps=2;
        
        % Prices, OAS, OAD & OAC
        prices=[];
        OAD=[];
        OAC=[];
        MOAD=[];
        MOAC=[];
        Carry=[];
        
        % CDF
        privateCDF=[];
        privateCDFgrid=(0:0.05:1)';
        
    end
    
    properties (Access = private, Dependent)
        % Annuity settings
        R_tilde=[];
        schedule=[];
        
        % FD time and space
        time=[];
        space=[];
        
        % FD settings
        J=[];
        dt=[];
        dr=[];
        
        % Distributions
        CDF=[];
    end
    
    properties (Access = public, Dependent)
        % Preallocate dependent model parameters
        kappa=[];
        sigma=[];
        theta=[];
        Xstar=[];
        R=[];
        
        % Dependent market
        market=[];
        
        % Time to maturity
        maturity=[];
        
        % FD settings
        Rannacher=true;
        smoothing=true;
        
        % OAS
        OAS=[];
    end
    
    methods (Access = public)
        function obj=cStantonHW(varargin)
            % Set user provided properties
            for i=1:2:nargin
                switch lower(varargin{i})
                    %%% Model parameters
                    case {'x','prepaymentcosts'}
                        % Set transaction costs
                        obj.X=varargin{i+1};
                    case {'f0','principal','notional'}
                        % Set principal
                        obj.F0=varargin{i+1};
                    case {'n','terms'}
                        % Set terms
                        obj.n=varargin{i+1};
                    case {'r','coupon'}
                        % Set coupon rate
                        obj.R=varargin{i+1};
                    case {'t','maturity'}
                        % Set maturity
                        obj.T=varargin{i+1};
                    case 'lambda'
                        % Set prepayment intensity
                        obj.lambda=varargin{i+1};
                    case {'a','alpha'}
                        % Set A
                        obj.Alpha=varargin{i+1};
                    case {'b','beta'}
                        % Set B
                        obj.Beta=varargin{i+1};
                    case 'market'
                        obj.mkt=varargin{i+1};
                        
                        % Finite difference parameters
                    case 'timesteps'
                        obj.n2=varargin{i+1};
                    case 'spacesteps'
                        obj.N=varargin{i+1};
                    case 'rannacher'
                        obj.Rannacher=varargin{i+1};
                    case 'smoothing'
                        obj.smoothing=varargin{i+1};
                    case 'feerateshort'
                        obj.feeRateShort=varargin{i+1};
                    case 'feeratelong'
                        obj.feeRateLong=varargin{i+1};
                end
            end
        end
        
        function [expectedReturn,accumulatedPrepayments,scenarios]=...
                expReturn(obj,shortRate,scenarios)
            % Get price today
            P0=obj.getPrice(shortRate);
            
            % Save maturity and Xstar
            saveT=obj.maturity;
            saveXstar=obj.Xstar;
            
            % Define scenarios
            if nargin<3
                scenarios=(-100:10:100)';
            end
            
            % Preallocate return vector
            expectedReturn=nan(size(scenarios));
            accumulatedPrepayments=nan(size(scenarios));
            
            % Loop over scenatios
            for j=1:size(scenarios,1)
                % Calculate X* for the four next draws and calc prices along
                % the way
                XstarForecast=nan(4,1);
                P=nan(4,1);
                for i=1:4
                    obj.maturity=saveT-0.25*i;
                    XstarForecast(i)=obj.getXstar(shortRate+scenarios(j)...
                        /10000);
                    obj.Xstar=[saveXstar;XstarForecast(1:i)];
                    P(i)=obj.getPrice(shortRate+scenarios(j)/10000);
                end
                
                % Prepayments
                temp=modelPrepayments(obj.Alpha,...
                    obj.Beta,...
                    obj.lambda,...
                    [saveXstar;XstarForecast],...
                    4);
                prepayments=temp(end-3:end);
                
                % Calculate capital evolution
                notional=100;
                for i=1:4
                    % Remaining payment dates
                    payDates=(saveT-i*0.25)*obj.n;
                    
                    % Calc payment
                    payment=obj.R_tilde/(1-(1+obj.R_tilde)^(-payDates))*...
                        notional;
                    
                    % Calc interest
                    interest=notional*obj.R_tilde;
                    
                    % Calc amortisation
                    amortisation=payment-interest;
                    
                    % Total draw
                    draw=amortisation*(1-prepayments(i))+...
                        notional*prepayments(i);
                    
                    % Total payment
                    totalPayment=interest+draw;
                    
                    % Subtract prepayments and amortisationfrom notional
                    notional=notional-draw;
                    
                    % Reinvest payment
                    notional=notional+totalPayment/P(i)*100;
                    
                end
                
                % Calculate return
                expectedReturn(j)=P(4)*notional/(P0*100)-1;
                
                % Calc accumulate prepayments
                accumulatedPrepayments(j)=1-prod(1-prepayments);
                
                % Reset maturity and Xstar
                obj.maturity=saveT;
                obj.Xstar=saveXstar;
            end   
        end
        
        function val=getPrice(obj,shortRate)
            % Check if prices are available
            if isempty(obj.prices)
                weightedPrices(obj);
            end
            
            % Hermite interpolate
            val=nan(size(shortRate));
            for i=1:size(val,1)
                val(i)=hermiteInterpolation(obj.space,obj.prices,...
                    shortRate(i));
            end
        end
        
        function [oad,oac,moad,moac,carry]=keyfigures(obj,shortRate)
            % Check if prices are available
            if isempty(obj.prices)
                weightedPrices(obj);
            end
            
            % Hermite interpolate keyfigures
            oad=nan(size(shortRate));
            for i=1:size(oad,1)
                oad(i)=hermiteInterpolation(obj.space,obj.OAD,...
                    shortRate(i));
            end
            if nargout>1
                oac=nan(size(shortRate));
                for i=1:size(oac,1)
                    oac(i)=hermiteInterpolation(obj.space,obj.OAC,...
                        shortRate(i));
                end
            end
            if nargout>2
                moad=nan(size(shortRate));
                for i=1:size(moad,1)
                    moad(i)=hermiteInterpolation(obj.space,obj.MOAD,...
                        shortRate(i));
                end
            end
            if nargout>3
                moac=nan(size(shortRate));
                for i=1:size(moac,1)
                    moac(i)=hermiteInterpolation(obj.space,obj.MOAC,...
                        shortRate(i));
                end
            end
            if nargout>4
                carry=nan(size(shortRate));
                for i=1:size(carry,1)
                    carry(i)=hermiteInterpolation(obj.space,obj.Carry,...
                        shortRate(i));
                end
            end
        end
        
        function val=getOAS(obj,shortRate,marketPrice)
            % Calc for OAS=0
            obj.OAS=0;
            weightedPrices(obj);
            modelPrice=hermiteInterpolation(obj.space,obj.prices,shortRate);
            
            
            % Set initial bounds
            lb=0;ub=0;
            if modelPrice>marketPrice
                while modelPrice>marketPrice
                    % Price at lower bound
                    Plb=modelPrice;
                    
                    % If upper bound is too low then increase it
                    ub=ub+0.05;
                    
                    % Calc for OAS=ub
                    obj.OAS=ub;
                    weightedPrices(obj);
                    modelPrice=hermiteInterpolation(obj.space,obj.prices,...
                        shortRate);
                    
                    % Price at upper bound
                    Pub=modelPrice;
                end
                lb=ub-0.05;
            else
                while modelPrice<marketPrice
                    % Price at upper bound
                    Pub=modelPrice;
                    
                    % If upper bound is too low then increase it
                    lb=lb-0.01;
                    
                    % Calc for OAS=ub
                    obj.OAS=lb;
                    weightedPrices(obj);
                    modelPrice=hermiteInterpolation(obj.space,obj.prices,...
                        shortRate);
                    
                    % Price at lb
                    Plb=modelPrice;
                end
                ub=lb+0.01;
            end
            
            deviation=1;
            while abs(deviation)>0.01
                %obj.OAS=0.5*(lb+ub);
                obj.OAS=lb+(ub-lb)/(Pub-Plb)*(marketPrice-Plb);
                weightedPrices(obj);
                modelPrice=hermiteInterpolation(obj.space,obj.prices,...
                    shortRate);
                deviation=modelPrice-marketPrice;
                if deviation>0
                    lb=obj.OAS;
                    Plb=modelPrice;
                else
                    ub=obj.OAS;
                    Pub=modelPrice;
                end
            end
            val=obj.OAS;
        end
        
        function val=getXstar(obj,shortRate,lb,ub)
            % if bounds are provided then use these
            if nargin<3
                lb=0;
                ub=1;
            end
            l=101;
            k=0;
            mb=0.5*(lb+ub);
            obj.X=mb;
            while abs(l-100*(1+obj.X))>0.01&&k<14
                mb=0.5*(lb+ub);
                obj.X=mb;
                [~,Ml,r]=obj.pricingFD;
                l=hermiteInterpolation(r,Ml,shortRate);
                if l>100*(1+obj.X)
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
            val=mb;
        end
        
        
        function [M_a,M_l,space,time]=pricingFD(obj)
            % Calculate annuity and payment dates
            annuitySchedule=obj.schedule;
            Y=annuitySchedule(:,4);
            fee=annuitySchedule(:,5)*obj.feeRateLong/obj.n;
            principal=annuitySchedule(:,5);
            
            % Terminal conditions
            M_a=zeros(obj.N+1,1);
            M_l=zeros(obj.N+1,1);
            
            % Get OAS
            if isempty(obj.OAS)
                spread=0;
            else
                spread=obj.OAS;
            end
            % Prepare for pricing
            obj.prepare;
            
            for j=obj.J+1:-1:1
                
                % Check if date is paydate
                if mod(j-1,obj.n2)==0
                    
                    % Get index for payment/principal
                    index=(j-1)/obj.n2+1;
                    
                    % OAS discounting
                    if j<(obj.J+1)
                        M_a=exp(-1/obj.n*spread)*M_a;
                    end
                    
                    % If paydate then calculate new principal
                    if j==1
                        F_t=obj.F0;
                    else
                        F_t=principal(index);
                    end
                    
                    % Delivery: min(F,M).  No-Delivery: F
                    if obj.delivery
                        G=min(F_t,M_a);
                    else
                        G=F_t;
                    end
                    
                    % Adjust for fees on 3M ARM
                    P=exp(-obj.space*obj.n);
                    G_l=(1+0.25*obj.feeRateShort*P).*(1+obj.X).*G;
                    
                    % Calculate Lambda (by smoothing)
                    b=G_l<M_l;b=b+0;
                    if obj.smoothing
                        input=M_l-G_l;
                        jumps=[0;diff(b)];
                        IX=(1:obj.N+1)'.*abs(jumps);IX(abs(jumps)==0)=[];
                        for i=1:size(IX,1)
                            minIX=max(IX(i)-25,1);
                            maxIX=min(IX(i)+24,size(b,1));
                            if minIX==1
                                b(minIX:maxIX)=...
                                    obj.smoothIndicator(...
                                    input(1:maxIX-minIX+1));
                            elseif maxIX==size(IX,1)
                                b(minIX:maxIX)=...
                                    obj.smoothIndicator(...
                                    input(size(b,1)-(maxIX-minIX)+1:end));
                            else
                                b(minIX:maxIX)=...
                                    obj.smoothIndicator(input(minIX:maxIX));
                            end
                        end
                    end
                    
                    LAMBDA=b.*(1-exp(-1/(obj.n)*obj.lambda));
                    
                    % Take "expected value"
                    if j>1
                        M_l=M_l.*(1-LAMBDA)+LAMBDA.*G_l;
                        M_a=M_a.*(1-LAMBDA)+LAMBDA.*G;
                    end
                    
                    % Add payment
                    M_a=M_a+Y(index);
                    M_l=M_l+Y(index)+fee(index);
                end
                
                % FD unless we are at t=0
                if j>1
                    % Calculate new values of liabilities
                    M_l=obj.matrixContainer{j,1}*M_l;
                    
                    % Calculate new values of bond
                    M_a=obj.matrixContainer{j,1}*M_a;
                end
            end
            
            % Set output
            space=obj.space;
            time=obj.time;
            
        end
        
        % Prepare for pricing.
        function prepare(obj)
            
            % Get forward curve and its derivative
            if isempty(obj.theta)
                forward=obj.mkt.forwardCurve(obj.time);
                derivative=obj.mkt.forwardCurveDiff(obj.time);
                
                % The derivative may explode due to strange behaviour of
                % the short end of the curve. Hence, we set the derivative
                % in [0;1) equal to the derivative at time t=1.
                B=obj.time<1;
                derivative(B)=derivative(sum(B));
                
                % Get theta
                obj.theta=forward+1/obj.kappa*derivative...
                    +obj.sigma^2/(2*obj.kappa)...
                    *(1-exp(-2*obj.kappa*obj.time));
            end
            
            % If FD matrices are empty then recalculate them
            if isempty(obj.matrixContainer)
                obj.updateFDmatrices;
            end
        end
        
    end
    
    % Private methods
    methods(Access = private)
        function weightedPrices(obj)
            
            % Get borrower distribution
            x=(obj.privateCDFgrid(2:end)+obj.privateCDFgrid(1:end-1))./2;
            probabilities=obj.CDF(2:end)-obj.CDF(1:end-1);
            
            % Get prices
            M=nan(obj.N+1,size(x,1));
            for i=1:size(x,1)
                obj.X=x(i);
                [M(:,i)]=pricingFD(obj);
            end
            % Weight prices
            obj.prices=M*probabilities;
            
            % Calc MOAD
            temp=-(obj.prices(3:end)-obj.prices(1:end-2))./(2*obj.dr*100);
            obj.OAD=[temp(1);temp;temp(end)];
            obj.MOAD=obj.OAD./obj.prices*100;
            
            % Calc MOAC
            temp=(obj.prices(3:end)-2*obj.prices(2:end-1)...
                +obj.prices(1:end-2))./(obj.dr^2*100^2);
            obj.OAC=[temp(1);temp;temp(end)];
            obj.MOAC=obj.OAC./obj.prices*100;
            
            % Get prices at t+1
            M=nan(obj.N+1,size(x,1));
            saveT=obj.T;obj.T=saveT-0.25;
            for i=1:size(x,1)
                obj.X=x(i);
                [M(:,i)]=pricingFD(obj);
            end
            obj.T=saveT;
            
            % Calc carry
            obj.Carry=(M*probabilities-obj.prices)./0.25;
            
        end
        
        function updateFDmatrices(obj)
            % Preallocation
            obj.matrixContainer=cell(obj.J+1,1);
            r=obj.space;
            
            for j=obj.J+1:-1:1
                if ((j-1)/obj.n2)==round((j-1)/obj.n2)&&obj.Rannacher
                    
                    % Perform Rannacher startup
                    dt05=obj.dt/obj.RannacherSteps;
                    
                    % Mu and sigma
                    Mu=obj.kappa*(obj.theta(j)-obj.space);
                    SIG=obj.sigma^2;
                    
                    % A B C D vectors
                    A=-1/(2*obj.dr)*Mu+1/2*SIG*1/(obj.dr^2);A=A(:);
                    B=-1/dt05-1/obj.dr^2*SIG-r;B=B(:);
                    C=1/(2*obj.dr)*Mu+1/(2*obj.dr^2)*SIG;C=C(:);
                    D=-1/dt05*ones(obj.N+1,1);D=D(:);
                    
                    % Correcting boundaries (zero convexity)
                    B(1)=-1/dt05-Mu(1)/obj.dr-r(1);
                    C(1)=Mu(1)/obj.dr;
                    D(1)=-1/dt05;
                    A(end)=-Mu(end)/obj.dr;
                    B(end)=-1/dt05+Mu(end)/obj.dr-r(end);
                    D(end)=-1/dt05;
                    
                    % Left Hand Side Matrix
                    LHS1=[[zeros(1,obj.N);diag(A(2:end))],zeros(obj.N+1,1)];
                    LHS2=diag(B);
                    LHS3=[zeros(obj.N+1,1),[diag(C(1:end-1));...
                        zeros(1,obj.N)]];
                    LHS=LHS1+LHS2+LHS3;
                    
                    % Right Hand Side Matrix
                    RHS=diag(D);
                    
                    % Invert LHS and multiply to RHS
                    obj.matrixContainer{j,1}=(LHS\RHS)^obj.RannacherSteps;
                else
                    % Mu and sigma
                    Mu=obj.kappa*(obj.theta(j)-obj.space);
                    SIG=obj.sigma^2;
                    
                    % A B C D vectors
                    A=1/(4*obj.dr)*Mu-1/(4*obj.dr^2)*SIG;A=A(:);
                    B=1/obj.dt+1/2*SIG*1/(obj.dr^2)+1/2*r;B=B(:);
                    C=-1/(4*obj.dr)*Mu-1/(4*obj.dr^2)*SIG;C=C(:);
                    D=1/obj.dt-1/(2*obj.dr^2)*SIG-1/2*r;D=D(:);
                    
                    % Correcting boundaries (zero convexity)
                    B(1)=1/obj.dt+Mu(1)/(2*obj.dr)+1/2*r(1);
                    C(1)=-Mu(1)/(2*obj.dr);
                    D(1)=1/obj.dt-Mu(1)/(2*obj.dr)-1/2*r(1);
                    A(end)=Mu(end)/(2*obj.dr);
                    B(end)=1/obj.dt-Mu(end)/(2*obj.dr)+1/2*r(end);
                    D(end)=1/obj.dt+Mu(end)/(2*obj.dr)-1/2*r(end);
                    
                    % Left Hand Side Matrix
                    LHS1=[[zeros(1,obj.N);diag(A(2:end))],zeros(obj.N+1,1)];
                    LHS2=diag(B);
                    LHS3=[zeros(obj.N+1,1),[diag(C(1:end-1));...
                        zeros(1,obj.N)]];
                    LHS=LHS1+LHS2+LHS3;
                    
                    % Right Hand Side Matrix
                    RHS1=[[zeros(1,obj.N);diag(-A(2:end))],...
                        zeros(obj.N+1,1)];
                    RHS2=diag(D);
                    RHS3=[zeros(obj.N+1,1),[diag(-C(1:end-1));...
                        zeros(1,obj.N)]];
                    
                    % Invert LHS and multiply to RHS
                    obj.matrixContainer{j,1}=LHS\(RHS1+RHS2+RHS3);
                end
            end
        end
        
        function y=smoothIndicator(~,x)
            one=ones(size(x));
            if all(x==0)
                y=one;
            else
                eps=(max(x)-min(x))/5;
                y=(0.5+1/(2*eps)*x+1/(2*pi)*...
                    sin(pi/eps*x)).*(x<=eps&-eps<=x)+one.*(x>eps);
            end
        end
        
        function clearFDmatrices(obj)
            obj.matrixContainer=[];
        end
        
        function clearPrices(obj)
            obj.prices=[];
            obj.MOAD=[];
            obj.MOAC=[];
        end
    end
    
    % Methods for private dependent properties.
    methods
        function val=get.R_tilde(obj)
            val=obj.R/obj.n;
        end
        
        function val=get.schedule(obj)
            val=annuity(obj.R,obj.F0,obj.T,obj.n);
        end
        
        function val=get.time(obj)
            val=(obj.t0:obj.dt:obj.T)';
        end
        
        function val=get.space(obj)
            val=(obj.r_min:obj.dr:obj.r_max)';
        end
        
        function val=get.J(obj)
            val=obj.T*obj.n*obj.n2;
        end
        
        function val=get.dt(obj)
            val=obj.T/obj.J;
        end
        
        function val=get.dr(obj)
            val=(obj.r_max-obj.r_min)/obj.N;
        end
        
        function set.CDF(obj,inputVal)
            obj.privateCDF=inputVal;
        end
        
        function val=get.CDF(obj)
            if ~isempty(obj.privateCDF)
                val=obj.privateCDF;
            else
                val=nan(size(obj.privateCDFgrid));
                for i=1:size(obj.privateCDFgrid,1)
                    [~,val(i)]=...
                        borrowerDistribution(obj.privateCDFgrid(i),...
                        obj.Xstar,obj.lambda,obj.Alpha,obj.Beta,obj.n);
                end
                obj.privateCDF=val;
            end
        end
    end
    
    % Get/Set for public properties
    methods
        function val=get.kappa(obj)
            val=obj.mkt.kappa;
        end
        function val=get.sigma(obj)
            val=obj.mkt.sigma;
        end
        function set.theta(obj,inputVal)
            obj.thetaPrivate=inputVal;
            obj.clearFDmatrices;
            obj.clearPrices;
        end
        function val=get.theta(obj)
            val=obj.thetaPrivate;
        end
        function set.OAS(obj,inputVal)
            obj.OASPrivate=inputVal;
            obj.clearPrices;
        end
        function val=get.OAS(obj)
            val=obj.OASPrivate;
        end
        function set.Xstar(obj,inputVal)
            obj.XstarPrivate=inputVal;
            obj.clearPrices;
            obj.privateCDF=[];
        end
        function val=get.Xstar(obj)
            val=obj.XstarPrivate;
        end
        function set.R(obj,inputVal)
            obj.OASPrivate=[];
            obj.RPrivate=inputVal;
            obj.clearPrices;
        end
        function val=get.R(obj)
            val=obj.RPrivate;
        end
        
        function set.market(obj,inputVal)
            obj.mkt=inputVal;
            obj.XstarPrivate=[];
            obj.thetaPrivate=[];
            obj.OASPrivate=[];
            obj.clearFDmatrices;
            obj.clearPrices;
        end
        
        function set.maturity(obj,inputVal)
            % Clear stuff if we select a longer maturity
            if inputVal>obj.T
                obj.thetaPrivate=[];
                obj.clearFDmatrices;
            end
            obj.clearPrices;
            obj.T=inputVal;
        end
        function val=get.maturity(obj)
            val=obj.T;
        end
        
        function set.Rannacher(obj,inputVal)
            obj.RannacherTimeStepping=inputVal;
            obj.clearFDmatrices;
            obj.clearPrices;
        end
        function val=get.Rannacher(obj)
            val=obj.RannacherTimeStepping;
        end
        
        function set.smoothing(obj,inputVal)
            obj.smoothPayoff=inputVal;
            obj.clearPrices;
        end
        function val=get.smoothing(obj)
            val=obj.smoothPayoff;
        end
    end
end

