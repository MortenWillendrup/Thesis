function createMarketData()
path='/users/helmig/dropbox/ku/speciale/Data/swap rates.xlsx';
[~,sheets]=xlsfinfo(path);

% Data structure
data=struct();

for i=1:size(sheets,2)
    % Get data
    a=xlsread(path,sheets{i});
    
    % Sort data
    [~,IX]=sort(a(:,1));
    a=a(IX,:);
    
    % Input to data structure
    if ~(strfind(sheets{i},'DKSW')==0)
        tenorString=sprintf('%sY',strrep(sheets{i},'DKSW',''));
        name=sprintf('SWAP%s',tenorString);
    elseif ~(strfind(sheets{i},'CIBO')==0)
        tenorString=strrep(sheets{i},'CIBO','');
        name=sprintf('CIBOR%s',tenorString);
    elseif ~(strfind(sheets{i},'EUCPAM')==0)
        tenorString=sprintf('%sY',strrep(sheets{i},'EUCPAM',''));
        name=sprintf('CAP%s',tenorString);
    elseif ~(strfind(sheets{i},'CITF')==0)
        tenorString=strrep(sheets{i},'CITF','');
        name=sprintf('CITA%s',tenorString);
    end
    
    % Convert tenor to double
    switch lower(right(tenorString,1))
        case 'w'
            tenor=str2double(left(tenorString,length(tenorString)-1))*1/52;
        case 'm'
            tenor=str2double(left(tenorString,length(tenorString)-1))*1/12;
        case 'y'
            tenor=str2double(left(tenorString,length(tenorString)-1));
    end
    
    % Save rates and dates
    dateVector=a(:,1)+datenum('1900-01-01')-2;
    data.(name).dates=(dateVector(1):dateVector(end))';
    [B,IX]=ismember(data.(name).dates,dateVector);
    data.(name).rates=nan(size(data.(name).dates));
    data.(name).rates(B)=a(IX(B),end);
    data.(name).tenorString=tenorString;
    data.(name).tenor=tenor;
end

% Save data
save('/users/helmig/dropbox/ku/speciale/Data/ratesContainer','data')
end