function [t, et]= solve_time_dependent_D(k, ek, initialT, eT, dT_dt, ttrials, D0, H)

if dT_dt==0 %isothermal
    A=k.^2/4;
    eA=k/2.*ek;
    D0 = DiffusionCoefficient(D0, H, Inf);
    D = DiffusionCoefficient(D0, H, initialT);
    t=A/D;
    ett=sqrt((eA./A).^2+(log(D/D0)*eT/initialT)^2)*100; %2 sigma percent uncertainty
    et=[-ett,ett]; %lower and upper
    return
end

% non-isothermal, use Monte Carlo method to estimate t
Tn=normrnd(initialT,eT,ttrials,1);
tn=NaN(ttrials,1);
t=NaN(length(k),1);
et=NaN(length(k),2);

MC_data=NaN(ttrials, length(k));
multiWaitbar( 'CloseAll' );
multiWaitbar(['Calculating for ', num2str(length(k)), ' k values'], 0, 'Color', [0.8 0.0 0.1]);

for i=1:length(k)
    multiWaitbar(['Calculating for ', num2str(length(k)), ' k values'], i/length(k));
    ki=k(i); eki=ek(i);
    kn=normrnd(ki,eki,ttrials,1);
    A=kn.^2/4;
    multiWaitbar('Monte Carlo model', 0, 'Color', [0.9 0.8 0.2]);
    for j=1:ttrials %Monte Carlo method
        multiWaitbar('Monte Carlo model', j/ttrials);
        syms a
        T=Tn(j)-dT_dt*a; %relationship between T and t
        D=matlabFunction(DiffusionCoefficient(D0, H, T));% or expressed directly as: D = @(t) (7e-8)*exp(-273000/8.314./(initialT-dT_dt*t));
        tn(j) = fzero(@(a)integral(D, 0, a)-A(j), 1);
    end
    MC_data(:,i)=tn;
    tn(isnan(tn))=[];
    logtn=log(tn);
    sd=std(logtn)*2;
    u=mean(logtn);
    t(i)=exp(u);
    et(i,1)=(exp(u-sd)/t(i)-1)*100;
    et(i,2)=(exp(u+sd)/t(i)-1)*100;
end
multiWaitbar( 'CloseAll' );

%histogram ln(t)
figure
sp = Subplots(length(k), -1);
for i=1:length(k)
    tp=MC_data(:,i);
    tp(isnan(tp))=[];
    logtp=log(tp);
    
    spi=sp.axis;
    histogram(logtp);
    xlabel('ln(t)')
    ylabel('Number')
    if isempty(logtp)
        spi.Title.String='No calculation result';
        continue
    end
    h=lillietest(logtp, 0.05);
    if h
        spi.Title.String=['Trial ', num2str(i), ', Total = ', num2str(length(logtp)), ', not normal distribution'];
    else
        spi.Title.String=['Trial ', num2str(i), ', Total = ', num2str(length(logtp)), ', normal distribution'];
    end
end

name1=repmat({'k='}, length(k), 1);
name2=mat2cell_wlg(k);
name=append(name1,name2)'; %or strcat
title=repmat({'t(s)'}, 1, length(k));

MC_output=num2cell(MC_data);
try
    dlmcell('Monte Carlo result.txt', name);
    dlmcell('Monte Carlo result.txt', title, '-a');
    dlmcell('Monte Carlo result.txt', MC_output, '-a');
catch
    filename=['Monte Carlo result_', datestr(now,'yyyy-mm-dd HH-MM-SS.FFF AM'), '.txt'];
    dlmcell(filename, name);
    dlmcell(filename, title, '-a');
    dlmcell(filename, MC_output, '-a');
end
