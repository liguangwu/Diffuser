function output = cell2num(celldata)
if ~iscell(celldata)
    errordlg('Please input cell array')
end

[r,c]=size(celldata);
output=NaN(r,c);
for i=1:r
    for j=1:c
        str=celldata{i,j};
        if isnumeric(str)
             output(i,j)=str;
        else
            output(i,j)=str2double(str);
        end
%         A=isstrprop(str,'digit');
%         B=str(A);
%         output(i,j)=str2double(B);
    end
end
end

