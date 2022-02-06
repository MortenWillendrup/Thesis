% ECBC figure

names={'Austria'
'Belgium'
'Denmark'
'Finland'
'France'
'Gernmany'
'Ireland'
'Italy'
'The Netherlands'
'Norway'
'Portugal'
'Spain'
'Sweden'
'Switzerland'
'United Kingdom'};

eurBn=[30894
16700
386232
33822
177813
207338
17062
138977
67604
113051
32970
232456
222444
117564
97127]./1000;


[eurBn,ix]=sort(eurBn,'descend');
names=names(ix);

f=figure('color','w');
b=bar(eurBn);
ax=gca;
ax.XAxis.TickLabels=names;
ax.XAxis.TickLabelRotation=35;
ax.XAxis.TickLabelInterpreter='latex';
ax.YAxis.TickLabelInterpreter='latex';
ax.YAxis.FontSize=12;
ax.XAxis.FontSize=12;