clc;clear;close all;
addpath('/Users/helmig/Dropbox/KU/Speciale/Data')

% Load market
mkt=market('date','2017-06-01');

mkt.swapCurve

[zero,forward,zeroTenors]=swapCalibration2(mkt.swapTenors,...
                                           mkt.swapCurve/100,...
                                           mkt.ciborTenors,...
                                           mkt.ciborCurve/100);
plot(mkt.swapTenors,mkt.swapCurve)
hold on
plot(zeroTenors,[zero,forward]*100)