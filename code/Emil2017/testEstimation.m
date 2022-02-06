% Clear
clc;clear;close all;

% Seed
seed=12345;
rng(seed)

% Time
T=120;

% Borrowers
N=10000;

% Some KNOWN function
t=(1:T)'./12;
Xstar=max(-0.0+cumsum(sqrt(0.001)*randn(size(t,1),1)),0);
%Xstar=0.4*sin(t*4)+0.45;
% Xstar=cumsum(ones(T,1)*0.05);
% Xstar=ones(T,1)*0.5;
% Parameters A and B to be estimated
A=2;
B=5;
X=betarnd(A,B,N,1);

% Prepayment intensity to be estimated per time dt
lambda=0.5;
dt=0.25;

% Prepayment indicator
tau=false(N,T);
U=rand(N,T);

% For each borrower (note: independent!)
for j=1:N
    % For each time unit
    for i=1:T
        % If this shit is fullfiled ...
        if X(j)<Xstar(i)
            % ... then prepayment happens with prob 1-exp(-lambda)
            if U(j,i)<1-exp(-lambda*dt)
                tau(j,i)=true;
                break;
            end
        end
    end
end

% Calculate prepayments
prepayments=sum(tau,1)';
Nt=N-cumsum(prepayments);
prepayRates=prepayments./Nt;

prepayProbs=nan(size(prepayRates));
for i=1:size(Xstar,1)
    prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,A,B,1/dt);
end


% Estimate
theta=[A;B;lambda];
lb=zeros(3,1);
ub=[10;10;2];
theta=fmincon(@(theta)GMMform2(theta,prepayRates,Xstar),theta,[],[],[],[],lb,ub);
A=theta(1);
B=theta(2);
lambda=theta(3);

% Extract expected prepayments
[~,expectedPrepay]=GMMform2(theta,prepayRates,Xstar);

% Plot prepayments
f=figure('color','w','position',[360   305   600   273]);
yyaxis left;
b=bar(t,prepayRates*100,...
      'barWidth',1,...
      'FaceColor',[0.8,0.8,0.8]);
hold on
plot(t,expectedPrepay*100,...
       'color',[0 0 0],...
       'linewidth',2)
ylabel('Prepayments (\%)','interpreter','latex')

% Plot costs
yyaxis right;
plot(t,Xstar*100,...
       'linestyle','--',...
       'color',[0 0 0])
ylabel('$X^*$ (\%)','interpreter','latex')
xlabel('Years','interpreter','latex')

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
xlim([0,10.1])

% Legend stuff
legend({'Prepayments (Left)\quad',...
        'Estimated Prepayment (Left)\quad',...
        '$X^*$ (Right)'},...
        'interpreter','latex',...
        'location','northoutside',...
        'orientation','horizontal',...
        'fontsize',12,...
        'box','off')

% Save
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\prepayGraph')

% Find distributions 
xGrid=(0:0.01:1)';
CDF=nan(size(xGrid));
PDF=nan(size(xGrid));
for i=1:size(xGrid,1)
    [PDF(i),CDF(i)]=borrowerDistribution(xGrid(i),Xstar,lambda,A,B,1/dt);
end

% Plot before and after
f=figure('color','w','position',[360   305   600   220]);
subplot(1,2,1)
% plot(xGrid,CDF)
% hold on
% [fx,x]=ecdf(X(sum(tau,2)==0));
% plot(x,fx)
yyaxis left
histogram(X,40,'FaceColor',[0.8,0.8,0.8])
ylabel('\# of borrowers','interpreter','latex')
yyaxis right
plot(xGrid,betapdf(xGrid,A,B),'color',[0 0 0],'linewidth',2)
ylabel('$f_t^X$','interpreter','latex')
xlabel('$X$','interpreter','latex')

% Layout 
ax=gca;
for i=1:2
    ax.YAxis(i).TickLabelInterpreter='latex';
    ax.YAxis(i).Label.Interpreter='latex';
    ax.YAxis(i).TickLength=[0,0];
    ax.YAxis(i).Label.FontSize=12;
    ax.YAxis(i).FontSize=12;
    ax.YAxis(i).Color=[0 0 0];
end
ax.XAxis.TickLabelInterpreter='latex';
ax.XAxis.Label.Interpreter='latex';
ax.XAxis.TickLength=[0,0];
ax.XAxis.Label.FontSize=10;
ax.XAxis.FontSize=10;
ax.XAxis.TickValues=[0 0.25 0.5 0.75 1];
ax.XAxis.TickLabels=[0 0.25 0.5 0.75 1]';
title('Initial distribution','interpreter','latex','fontsize',12)

subplot(1,2,2)
% plot(xGrid,CDF)
% hold on
% [fx,x]=ecdf(X(sum(tau,2)==0));
% plot(x,fx)
yyaxis left
histogram(X(sum(tau,2)==0),100,'FaceColor',[0.8,0.8,0.8])
ylabel('\# of borrowers','interpreter','latex')
yyaxis right
plot(xGrid,PDF,'color',[0 0 0],'linewidth',2)
ylabel('$f_t^X$','interpreter','latex')
xlabel('$X$','interpreter','latex')

% Layout 
ax=gca;
for i=1:2
    ax.YAxis(i).TickLabelInterpreter='latex';
    ax.YAxis(i).Label.Interpreter='latex';
    ax.YAxis(i).TickLength=[0,0];
    ax.YAxis(i).Label.FontSize=12;
    ax.YAxis(i).FontSize=12;
    ax.YAxis(i).Color=[0 0 0];
end
ax.XAxis.TickLabelInterpreter='latex';
ax.XAxis.Label.Interpreter='latex';
ax.XAxis.TickLength=[0,0];
ax.XAxis.Label.FontSize=10;
ax.XAxis.FontSize=10;
ax.XAxis.TickValues=[0 0.25 0.5 0.75 1];
ax.XAxis.TickLabels=[0 0.25 0.5 0.75 1]';
title('Final distribution','interpreter','latex','fontsize',12)

% Save
%saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\prepayDist')



function [val,f]=GMMform(theta,prepayments,Xstar,N)
% Get parameters
A=theta(1);
B=theta(2);
lambda=theta(3);
dt=0.25;

% Calculate Beta probabilities
BETA=betacdf(Xstar,A,B);

% Calculate survival probabilities
survival=cumprod([1;(1-(1-exp(-lambda*dt))*BETA(2:end))]);

% Calculate prepayment probabilities
P=(1-exp(-lambda*dt))*BETA.*survival;

% Theoretical moments
f=P.*N;
%f=P./(1-cumsum([0;P(1:end-1)]));

% Calculate quadratic form
val=(f-prepayments)'*(f-prepayments);
end



function [val,expectedPrepayRates]=GMMform2(theta,prepayments,Xstar)
% Get parameters
A=theta(1);
B=theta(2);
lambda=theta(3);
dt=0.25;

% Calculate prepayment probabilities
prepayProbs=nan(size(prepayments));
for i=1:size(Xstar,1)
    prepayProbs(i)=prepaymentProbability(Xstar(i),Xstar(1:i-1),lambda,A,B,1/dt);
end

% Calculate prepayment rates
remaining=1-cumsum(prepayProbs);
expectedPrepayRates=prepayProbs./remaining;

% Calculate quadratic form
f=prepayments-expectedPrepayRates;
val=f'*f;
end