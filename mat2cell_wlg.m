function cdata = mat2cell_wlg(data)
%将数值矩阵转化为元胞，如果是NaN或Inf，转化为''
[r,c]=size(data);
cdata=cell(r,c);
for i=1:r
    for j=1:c
        if isnan(data(i,j)) || isinf(data(i,j))
            cdata{i,j}='';
        else
            cdata{i,j}=num2str(data(i,j));
%             cdata{i,j}=data(i,j);
        end
    end
end
end

