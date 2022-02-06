classdef cStanton
    % Pricing in the Stanton (1995) model
    
    properties (Access = public)
        % Preallocate model parameters
        kappa=[];
        mu=[];
        sig=[];
        q=[];
        X=[];
        F0=[];
        R=[];
        n=[];
        T=[];
        lambda1=[];
        lambda2=[];
        
        % Preallocate finite difference settings
        % Set discretisation parameters
        n2=4;
        J=[];
        dt=[];
        t0=0;
        r_min=0.0;
        r_max=0.9;
        N=400;
        dr=[];
    end
    
    properties(Dependent)
        R_tilde=[];
        schedule=[];
    end
    
    methods (Access = public)
        function obj=cStanton(varargin)
            % Set user provided properties
            for i=1:2:nargin
                switch lower(varargin{i})
                    %%% Model parameters
                    case 'kappa'
                        % Set kappa
                        obj.kappa=varargin{i+1};
                    case 'mu'
                        % Set mu
                        obj.mu=varargin{i+1};
                    case {'sig','sigma'}
                        % Set sigma
                        obj.sig=varargin{i+1};
                    case 'q'
                        % Set q (market price of risk)
                        obj.q=varargin{i+1};
                    case 'x'
                        % Set transaction costs
                        obj.X=varargin{i+1};
                    case {'f0','principal'}
                        % Set principal
                        obj.F0=varargin{i+1};
                    case {'n','terms'}
                        % Set principal
                        obj.n=varargin{i+1};
                    case {'r','coupon'}
                        % Set coupon rate
                        obj.R=varargin{i+1};
                    case {'t','maturity'}
                        % Set maturity
                        obj.T=varargin{i+1};
                    case 'lambda1'
                        % Set maturity
                        obj.lambda1=varargin{i+1};
                    case 'lambda2'
                        % Set maturity
                        obj.lambda2=varargin{i+1};
                        
                    % Finite difference parameters
                    
                end
            end

            % Calculate finite difference settings
            % Set discretisation parameters
            obj.J=obj.T*obj.n*obj.n2;
            obj.dt=obj.T/obj.J;
            obj.dr=(obj.r_max-obj.r_min)/obj.N;
        end
        
        function [M_a,M_l,space,time]=pricingFD(obj)
            % Calculate annuity and payment dates
            annuitySchedule=obj.schedule;
            Y=annuitySchedule(:,4);
            principal=annuitySchedule(:,5);
            
            % Create grid
            time=(obj.t0:obj.dt:obj.T)';
            space=(obj.r_min:obj.dr:obj.r_max)';
            
            % Terminal conditions
            M_a=zeros(obj.N+1,1);
            M_l=zeros(obj.N+1,1);

            % Mu and Sigma
            Mu=obj.kappa*obj.mu-(obj.kappa+obj.q)*space;
            Sig=obj.sig^2*space;
            
            % A B C D
            A=1/(4*obj.dr)*Mu-1/(4*obj.dr^2)*Sig;A=A(:);
            B=1/obj.dt+1/2*Sig*1/(obj.dr^2)+1/2*space;B=B(:);
            C=-1/(4*obj.dr)*Mu-1/(4*obj.dr^2)*Sig;C=C(:);
            D=1/obj.dt-1/(2*obj.dr^2)*Sig-1/2*space;D=D(:);

            % Correcting boundaries (zero convexity)
            B(1)=1/obj.dt+Mu(1)/(2*obj.dr)+1/2*space(1);
            C(1)=-Mu(1)/(2*obj.dr);
            D(1)=1/obj.dt-Mu(1)/(2*obj.dr)-1/2*space(1);
            A(end)=Mu(end)/(2*obj.dr);
            B(end)=1/obj.dt-Mu(end)/(2*obj.dr)+1/2*space(end);
            D(end)=1/obj.dt+Mu(end)/(2*obj.dr)-1/2*space(end);

            % Left Hand Side Matrix
            LHS1=[[zeros(1,obj.N);diag(A(2:end))],zeros(obj.N+1,1)];
            LHS2=diag(B);
            LHS3=[zeros(obj.N+1,1),[diag(C(1:end-1));zeros(1,obj.N)]];
            LHS=LHS1+LHS2+LHS3;

            % Right Hand Side Matrix
            RHS1=[[zeros(1,obj.N);diag(-A(2:end))],zeros(obj.N+1,1)];
            RHS2=diag(D);
            RHS3=[zeros(obj.N+1,1),[diag(-C(1:end-1));zeros(1,obj.N)]];
            RHS=RHS1+RHS2+RHS3;

            % Invert LHS
            invLHS=inv(LHS);

            for j=obj.J+1:-1:1
                
                % Check if date is paydate
                if mod(j-1,obj.n2)==0
                                        
                    % Get index for payment/principal
                    index=(j-1)/obj.n2+1;
                    
                    % If paydate then calculate new principal
                    if j==1
                        F_t=obj.F0;
                    else
                        F_t=principal(index);
                    end
                    
                    % Calculate Lambda
                    b=M_l>(1+obj.X)*F_t;
                    LAMBDA=obj.lambda1+obj.lambda2*b;
                    LAMBDA2=1-exp(-1/(obj.n)*LAMBDA);

                    % Take "expected value"
                    M_l=M_l.*(1-LAMBDA2)+LAMBDA2.*F_t*(1+obj.X);
                    M_a=M_a.*(1-LAMBDA2)+LAMBDA2.*F_t;

                    % Add payment
                    M_a=M_a+Y(index);
                    M_l=M_l+Y(index);
                end
                
                % Calculate new values of liabilities
                f=RHS*M_l;
                M_l=(invLHS)*f;

                % Calculate new values of bond
                f=RHS*M_a;
                M_a=(invLHS)*f;
            end            
        end
    end
    
    % Methods for dependent properties. 
    % All dynamic proterties have empty set functions
    methods 
        function val=get.R_tilde(obj)
            val=obj.R/obj.n;
        end
        function obj=set.R_tilde(obj,~);end;
        
        function val=get.schedule(obj)
            val=annuity(obj.R,obj.F0,obj.T,obj.n);
        end
        function obj=set.schedule(obj,~);end;
    end
    
end

