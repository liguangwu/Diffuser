function [name, D0, H] = readDiffCoef()

str='DiffusionCoefficient.xlsx';
if exist(str,'file')
    [~, sheetname] = xlsfinfo(str);
    [~,cs]=size(sheetname);
    if cs==1
        sheet=sheetname{1};
        v=1;
    else
        [st,v]=listdlg('PromptString','Select a sheet','SelectionMode','single','ListString',sheetname);
        sheet=sheetname{st};
    end
    if v
        [~, ~, raw]=xlsread(str,sheet);
    else
        name=[]; D0=[]; H=[];
        return
    end
else
    errordlg("'DiffusionCoefficient.xlsx' does not exist, please select a diffusion coefficient excel manually.")
    name=[]; D0=[]; H=[];
    return
end
name=join(raw(2:end,1:2));
D0=raw(2:end,3);
D0=cell2mat(D0);
H=raw(2:end,4);
H=cell2mat(H);