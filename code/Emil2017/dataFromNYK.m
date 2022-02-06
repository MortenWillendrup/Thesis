isin='DK0009753469';
url=sprintf('https://www.nykredit.com/investorinfo/GetFileServlet?filename=%s.xls',isin);
dataPath=urlwrite(url,'/Users/helmig/Dropbox/KU/Speciale/Kode/tesdata.xls');
[~,sheets]=xlsfinfo(dataPath);
[a,b,c]=xlsread(dataPath,sheets{4},'A5:F1000');
