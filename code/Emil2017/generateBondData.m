pathXLS='/users/helmig/dropbox/ku/speciale/Data/prepayments.xls';
[~,sheets]=xlsfinfo(pathXLS);
[~,~,info]=xlsread(pathXLS,'Info');

bonds=struct();
bonds.isins=sheets(2:end)';
for i=2:size(sheets,2)
    bonds.(sheets{i}).Bank=info{i,2};
    bonds.(sheets{i}).Coupon=info{i,3};
    bonds.(sheets{i}).Maturity=info{i,4}-49583+datenum('2035-10-01');
    [~,~,temp]=xlsread(pathXLS,sheets{i});temp(1,:)=[];
    [dates,IX]=sort(datenum(temp(:,1),'dd-mm-yyyy'));
    bonds.(sheets{i}).Dates=dates;
    bonds.(sheets{i}).Ordinary=cell2mat(temp(IX,2));
    bonds.(sheets{i}).Extraordinary=cell2mat(temp(IX,3));
    bonds.(sheets{i}).TotalDraw=cell2mat(temp(IX,4));
    bonds.(sheets{i}).Outstanding=cell2mat(temp(IX,5));
end

save('/users/helmig/dropbox/ku/speciale/Data/bonds','bonds')