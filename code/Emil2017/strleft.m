function strOut=strleft(string,numChar)
strOut=cell(size(string));
string=string(:);
for i=1:size(string,1)
    strOut{i}=left(string{i},numChar);
end
end