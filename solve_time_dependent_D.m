function [t, et, isbreak]= solve_time_dependent_D(k, ek, initialT, eT, dT_dt, ttrials, DiffCoef, coolingpath, Dtype, app)
%ek, and eT are 2sigma absolute errors, et is 2 sigma percent error
isbreak=0;
if Dtype %use Arrhenius equation
    D0=DiffCoef{1};
    H=DiffCoef{2};
else %use a constant D
    D=DiffCoef(1);
    eD=DiffCoef(2);
    if isnan(eD); eD=0; end
end

if dT_dt==0 || ~Dtype%isothermal
    Dt=10.^k;
    eDt=log(10)*(10.^k).*ek;
    if Dtype
       D = DiffusionCoefficient(D0, H, initialT);
       ett=sqrt((eDt./Dt).^2+(log(D/D0)*eT/initialT)^2)*100; %2 sigma percent error
    else
       ett=sqrt((eDt./Dt).^2+(eD/D)^2)*100; %2 sigma percent error
    end
    t=Dt/D;
    et=[-ett,ett]; %lower and upper 2 sigma percent error
    return
end

% non-isothermal, use Monte Carlo method to estimate t
Tn=normrnd(initialT,eT,ttrials,1);
tn=NaN(ttrials,1);
t=NaN(length(k),1);
et=NaN(length(k),2);

MC_data=NaN(ttrials, length(k));

wb=uiprogressdlg(app.DiffuserGUI,'ShowPercentage', 'on', 'Cancelable', 'on');

for i=1:length(k)
    if isbreak
        break
    end
    wb.Title=['Calculating for ', num2str(length(i)), '/', num2str(length(k)), ' Dt value'];
    ki=k(i); eki=ek(i);
    kn=normrnd(ki,eki,ttrials,1);
    Dt=10.^kn;
    for j=1:ttrials %Monte Carlo method
        if wb.CancelRequested
            isbreak=1;
            break
        end
        wb.Message='Monte Carlo estimating time uncertainty';
        wb.Value=j/ttrials;
        syms a
        %relationship between T and t
        switch coolingpath
            case 'isothermal'
                T=Tn(j);
            case 'linear cooling'
                T=Tn(j)-dT_dt*a;
            case 'exponential cooling'
                T=Tn(j)*exp(-dT_dt*a);
            case 'parabolic cooling'
                T=Tn(j)-dT_dt*(a^2);
        end
        D=matlabFunction(DiffusionCoefficient(D0, H, T));% or expressed directly as: D = @(t) (7e-8)*exp(-273000/8.314./(initialT-dT_dt*t));
        tn(j) = fzero(@(a)integral(D, 0, a)-Dt(j), 1);
    end
    MC_data(:,i)=tn;
    tn(isnan(tn))=[];
    logtn=log(tn);
    sd=std(logtn)*2; %2 sigma absoute error
    u=mean(logtn);
    t(i)=exp(u);
    et(i,1)=(exp(u-sd)/t(i)-1)*100; %2 sigma percent error
    et(i,2)=(exp(u+sd)/t(i)-1)*100;
end
close(wb);
if isbreak
    return
end

%histogram ln(t)
normaltest=cell(1,length(k));
figure('Name', 'ln(t) distribution when cooling')
sp = Subplots(length(k), -1);
for i=1:length(k)
    tp=MC_data(:,i);
    tp(isnan(tp))=[];
    logtp=log(tp);
    sd=std(logtp)*2; %2 sigma absoute error
    u=mean(logtp);
    spi=sp.axis;
    histogram(logtp);
    xlabel('ln(t)')
    ylabel('Number')
    if isempty(logtp)
        spi.Title.String='No calculation result';
        normaltest(i)={'not normal distribution'};
        continue
    end
    h=lillietest(logtp, 0.05);
    if h
        legend(['log10(Dt) = ', num2str(k(i),'%.2f'), newline, ...
            'ln(t) = ', num2str(u,'%.2f'), ' ± ', num2str(sd,'%.2f'), ' (2SD, N = ', num2str(length(logtp)), ')', newline, ...
            'not normal distribution']);
        normaltest(i)={'not normal distribution'};
    else
        legend(['log10(Dt) = ', num2str(k(i),'%.2f'), newline, ...
            'ln(t) = ', num2str(u,'%.2f'), ' ± ', num2str(sd,'%.2f'), ' (2SD, N = ', num2str(length(logtp)), ')', newline, ...
            'normal distribution']);
        normaltest(i)={'normal distribution'};
    end
end

name1=repmat({'log10(Dt)='}, length(k), 1);
name2=cell(length(k),1);
for i=1:length(k)
        name2{i}=num2str(k(i));
end
name=append(name1,name2)'; %or strcat
str1=repmat({'t(s), '}, 1, length(k));
tstr=append(str1, normaltest);
MC_output=num2cell(MC_data);

%write into file
filename=['Monte Carlo result ',datestr(now,'yyyy-mm-dd HH-MM-SS.FFF AM'), '.txt'];
[FileName,PathName,~] = uiputfile({'*.txt';'*.csv';'*.xls';'*.xlsx'},'Save result',filename);
if ~FileName
    return
end
[~,~,extension]=fileparts(FileName);
switch extension
    case {'.txt','.csv'}
        writecell([name;tstr;MC_output],fullfile(PathName,FileName),"Delimiter",'tab');
    case {'.xlsx','.xls'}
        %                                         writecell([name;tstr;MC_output],fullfile(PathName,FileName),'AutoFitWidth', 0);
        xlswrite(fullfile(PathName,FileName),[name;tstr;MC_output]);
    otherwise
        return
end