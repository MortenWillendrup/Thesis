%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   Prepayment figure
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;

rho=[0.5;1;3;5;10;1000];
dt=1;
time=(0.00:1/12:5)';
pool=nan(size(time));
pool(1)=1;
RHO=[0.5,1,3,5,10,1000];

% Delta is a CONSTANT fraction of the pool for whom it is optimal to prepay
delta=0.5;
%
for j=1:size(RHO,2)
for i=1:size(time,1)
    rho=RHO(j);
    remaining=exp(-rho*delta*time(i));
    pool(i)=(1-exp(-rho*delta))*remaining;
end
hold on
plot(time*12,pool)
end



