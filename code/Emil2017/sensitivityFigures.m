%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Sensitivity figures
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

% Load market
endDate=datenum('2017-06-01');
mkt=market('date',endDate);

% Parameters
F0=100;
R=0.025;
n=4;
T=30;
lambda=0.5169;
Alpha=0.1107;
Beta=2.3877;

% Initiate model
model=cStantonHW('F0',F0,...
                 'R',R,...
                 'n',n,...
                 'T',T,...
                 'lambda',lambda,...
                 'A',Alpha,...
                 'B',Beta,...
                 'market',mkt,...
                 'timesteps',3,...
                 'spacesteps',400,...
                 'rannacher',false,...
                 'smoothing',false);
             
% Regular Crank-Nicolson 
rGrid=(-0.05:0.00001:0.1)';
model.T=30;
[OAD,OAC]=model.keyfigures(rGrid);

% Draw a figure
f=figure('color','w','position',[360   305   600   230]);
s=subplot(1,2,1);
latexPlot('f',s,'x',rGrid*100,'y',OAD,...
          'xlabel','Short Rate (\%)',...
          'ylabel','OAD')
s.Position(2)=0.17;s.Position(4)=0.75;
s=subplot(1,2,2);
latexPlot('f',s,'x',rGrid*100,'y',OAC,...
          'xlabel','Short Rate (\%)',...
          'ylabel','OAC')
s.Position(2)=0.17;s.Position(4)=0.75;
saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\OADOAC')

% Rannacher time stepping
model.Rannacher=true;
[OAD,OAC]=model.keyfigures(rGrid);

% Draw a figure
f=figure('color','w','position',[360   305   600   230]);
s=subplot(1,2,1);
latexPlot('f',s,'x',rGrid*100,'y',OAD,...
          'xlabel','Short Rate (\%)',...
          'ylabel','OAD')
s.Position(2)=0.17;s.Position(4)=0.75;
s=subplot(1,2,2);
latexPlot('f',s,'x',rGrid*100,'y',OAC,...
          'xlabel','Short Rate (\%)',...
          'ylabel','OAC')
s.Position(2)=0.17;s.Position(4)=0.75;
saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\OADOAC2')

% Rannacher time stepping + smoothing
model.smoothing=true;
[OAD,OAC]=model.keyfigures(rGrid);

% Draw a figure
f=figure('color','w','position',[360   305   600   230]);
s=subplot(1,2,1);
latexPlot('f',s,'x',rGrid*100,'y',OAD,...
          'xlabel','Short Rate (\%)',...
          'ylabel','OAD')
s.Position(2)=0.17;s.Position(4)=0.75;
s=subplot(1,2,2);
latexPlot('f',s,'x',rGrid*100,'y',OAC,...
          'xlabel','Short Rate (\%)',...
          'ylabel','OAC')
s.Position(2)=0.17;s.Position(4)=0.75;
hold on
% Fill
% XXX=[min(rGrid(OAC<0));max(rGrid(OAC<0))];XXX=[XXX;flip(XXX)];
% YYY=[repmat(-2,size(XXX,1)/2,1);repmat(1+0,size(XXX,1)/2,1)];
% ff=fill(XXX*100,YYY,[0 0 0]);
% ff.FaceAlpha=0.1;
% ff.EdgeAlpha=0;
saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\OADOAC3')