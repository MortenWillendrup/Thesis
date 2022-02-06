%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Beta dist
%   
%   Emil Rode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
x=(0:0.00001:1)';

Alpha=[1 3 2 1.5]';
Beta=[3 1 2 5]';

y=nan(size(x,1),size(Alpha,1));
l=cell(4,1);
for i=1:size(Alpha,1)
    y(:,i)=betapdf(x,Alpha(i),Beta(i));
    l{i}=sprintf('$\\alpha=%s, \\beta=%s$\quad',num2str(Alpha(i),'%4.1f'),...
        num2str(Beta(i),'%4.1f'));
end


f=latexPlot('x',x,'y',y,'legend',l,'location','eastoutside','orientation','vertical')
ylim([0 3])
saveFigAsPdf(f,'\users\helmig\dropbox\KU\speciale\thesis\betaDist')