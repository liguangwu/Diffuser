function [t, et, MC_data, error_contribution, isbreak]= solve_time_dependent_D(k, ek, initialT, eT, dT_dt, ttrials, DiffCoef, coolingpath, Dtype, withD0Eerror, app)
%ek, and eT are 2sigma absolute errors, et is 2 sigma percent error [-2s%, +2s%]
isbreak=0; MC_data=[]; error_contribution=[];
%use a constant D-------------------------------
if ~Dtype
    D=DiffCoef(1);
    eD=DiffCoef(2); %2 sigma
    if isnan(eD); eD=0; end
    Dt=10.^k;
    eDt=log(10)*(10.^k).*ek;
    elnt=sqrt((eDt./Dt).^2+(eD/D)^2); %2 sigma, σln(t)=σt/t
    t=Dt/D;
    lnt=log(t);
    et(:,1)=(exp(lnt-elnt)./t-1)*100; %2 sigma percent error
    et(:,2)=(exp(lnt+elnt)./t-1)*100;
    return
end
%------------------------------------------------
%use Arrhenius equation, use Monte Carlo method to estimate t
lnD0=DiffCoef{1}; D0=exp(lnD0);
elnD0=DiffCoef{2}/2; %1 sigma
if isnan(elnD0); elnD0=0; end
E=DiffCoef{3}; %kJ/mol
eE=DiffCoef{4}/2;  %1 sigma
if isnan(eE); eE=0; end
Cov=DiffCoef{5}; %covariance between ln(D0) and E
if isnan(Cov); Cov=0; end
%1 = error not considered, 2 = error considered
cov1=[0 0; 0 0];
cov2=[elnD0^2 Cov; Cov eE^2]; %covariance matrix between ln(D) and E
covmatrix1=[cov1 [0;0]; [0,0] 0]; %covariance matrix between ln(D), E and T
covmatrix2=[cov1 [0;0]; [0,0] eT^2/4]; %covariance matrix between ln(D), E and T
covmatrix3=[cov2 [0;0]; [0,0] eT^2/4]; %covariance matrix between ln(D), E and T
Tn1=normrnd(initialT,0,ttrials,1);
Tn2=normrnd(initialT,eT/2,ttrials,1);
R1 = mvnrnd([lnD0,E],cov1,ttrials);
D0n1=exp(R1(:,1)); En1=R1(:,2); 
R2 = mvnrnd([lnD0,E],cov2,ttrials);
D0n2=exp(R2(:,1)); En2=R2(:,2);
t=NaN(length(k),1);
et=NaN(length(k),2);
MC_data=NaN(ttrials, length(k)*2);
% isothermal------------------------------------
if dT_dt==0
    Dt=10.^k;
    eDt=log(10)*(10.^k).*ek; %2 sigma
    D=D0*exp(-E*1000/8.314/initialT);
    t=Dt/D;
    lnt=log(t);
    J=D*[1;-1000/8.314/initialT;E*1000/8.314/(initialT^2)]; %Jacobi matrix between ln(D), E and T
    eD1=sqrt(J'*covmatrix1*J)*2; %2 sigma, %only consider curve-fitting error
    eD2=sqrt(J'*covmatrix2*J)*2; %2 sigma, %consider curve-fitting and temperature errors
    eD3=sqrt(J'*covmatrix3*J)*2; %2 sigma, %consider curve-fitting, temperature, and experimental errors
    elnt1=sqrt((eDt./Dt).^2+(eD1/D)^2); %2 sigma, σln(t)=σt/t
    elnt2=sqrt((eDt./Dt).^2+(eD2/D)^2); %2 sigma, σln(t)=σt/t
    elnt3=sqrt((eDt./Dt).^2+(eD3/D)^2); %2 sigma, σln(t)=σt/t
    e1=elnt1./elnt3*100; e2=(elnt2./elnt3-elnt1/elnt3)*100; e3=(1-elnt2./elnt3)*100; %contribution to ln(t) error from curve fitting, temperature, and experimental errors
    if withD0Eerror
        et(:,1)=(exp(lnt-elnt3)./t-1)*100; %2 sigma percent error
        et(:,2)=(exp(lnt+elnt3)./t-1)*100;
        error_contribution=[e1(1,1),e2(1,1),e3(1,1)];
    else
        et(:,1)=(exp(lnt-elnt2)./t-1)*100; %2 sigma percent error
        et(:,2)=(exp(lnt+elnt2)./t-1)*100;
        e11=e1./(e1+e2)*100;
        e22=e2./(e1+e2)*100;
        error_contribution=[e11(1,1),e22(1,1),0];
    end
    %generate random T and corresponding t
    for i=1:length(k)
        Tn2=normrnd(initialT,eT/2,ttrials,1);
        if withD0Eerror
            Dn = D0n2.*exp(-En2*1000/8.314./Tn2);
        else
            Dn = D0n1.*exp(-En1*1000/8.314./Tn2);
        end
        kn=normrnd(k(i),ek(i)/2,ttrials,1);
        Dtn=10.^kn;
        MC_data(:,i*2-1)=Tn2-273.15;
        MC_data(:,i*2)=Dtn./Dn;
    end
end
% ==============================================================
%non-isothermal------------------------------------------------ 
if dT_dt~=0
    tn=NaN(ttrials,1); tn1=NaN(ttrials,1); tn2=NaN(ttrials,1); tn3=NaN(ttrials,1);
    wb=uiprogressdlg(app.DiffuserGUI,'ShowPercentage', 'on', 'Cancelable', 'on');
    if length(k)==1
        i=1;
        MC_data(:,i*2-1)=Tn2-273.15;
        wb.Title=['Calculating for ', num2str(i), '/', num2str(length(k)), ' Dt value'];
        ki=k(i); eki=ek(i);
        kn=normrnd(ki,eki/2,ttrials,1);
        Dt=10.^kn;
        for j=1:ttrials %Monte Carlo method
            if wb.CancelRequested
                isbreak=1;
                break
            end
            wb.Message='Monte Carlo estimating time uncertainty';
            wb.Value=j/ttrials;
            %relationship between T and t
            switch coolingpath
                case 'isothermal'
                    D1 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-0*a)); %only consider curve-fitting error
                    D2 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-0*a)); %consider curve-fitting and temperature errors
                    D3 =@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-0*a)); %consider curve-fitting, temperature, and experimental errors
                case 'linear cooling'
                    D1 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-dT_dt*a));
                    D2 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-dT_dt*a));
                    D3 =@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-dT_dt*a));
                case 'exponential cooling'
                    D1 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)*exp(-dT_dt*a)));
                    D2 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)*exp(-dT_dt*a)));
                    D3 =@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)*exp(-dT_dt*a)));
                case 'parabolic cooling'
                    D1 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-dT_dt*(a^2)));
                    D2 =@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-dT_dt*(a^2)));
                    D3 =@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-dT_dt*(a^2)));
            end
            tn1(j) = fzero(@(a)integral(D1, 0, a)-Dt(j), 1);
            tn2(j) = fzero(@(a)integral(D2, 0, a)-Dt(j), 1);
            tn3(j) = fzero(@(a)integral(D3, 0, a)-Dt(j), 1);
        end
        if withD0Eerror
            MC_data(:,i*2)=tn3;
            tn1(isnan(tn1))=[]; tn2(isnan(tn2))=[]; tn3(isnan(tn3))=[];
            logtn1=log(tn1); logtn2=log(tn2); logtn3=log(tn3);
            sd1=std(logtn1)*2; sd2=std(logtn2)*2; sd3=std(logtn3)*2; %2 sigma absoute error
            u1=mean(logtn1); u2=mean(logtn2); u3=mean(logtn3);
            t(i)=exp(u3);
            et(i,1)=(exp(u3-sd3)/t(i)-1)*100; %2 sigma percent error
            et(i,2)=(exp(u3+sd3)/t(i)-1)*100;
        else
            MC_data(:,i*2)=tn2;
            tn1(isnan(tn1))=[]; tn2(isnan(tn2))=[];
            logtn1=log(tn1); logtn2=log(tn2);
            sd1=std(logtn1)*2; sd2=std(logtn2)*2; %2 sigma absoute error
            u1=mean(logtn1); u2=mean(logtn2);
            t(i)=exp(u2);
            et(i,1)=(exp(u2-sd2)/t(i)-1)*100; %2 sigma percent error
            et(i,2)=(exp(u2+sd2)/t(i)-1)*100;
            u3=u2; sd3=sd2;
        end
        sd1=sd1+u3-u1; sd2=sd2+u3-u2; %keep the average same
        if sd1<=sd2 && sd2<=sd3 %otherwise, the trials are not enough
            e1=sd1/sd3*100; e2=(sd2/sd3-sd1/sd3)*100; e3=(1-sd2/sd3)*100; %contribution to ln(t) error from curve fitting, temperature, and experimental errors
            error_contribution=[e1,e2,e3];
        end
    end
    %------------------------------------------------ 
    %calculate t when no flat peak/trough in profiles I-L
    if length(k)>1
        if withD0Eerror
            D0n=D0n2; En=En2;
        else
            D0n=D0n1; En=En1;
        end
        for i=1:length(k)
            if isbreak
                break
            end
            Tn=normrnd(initialT,eT/2,ttrials,1);
            MC_data(:,i*2-1)=Tn-273.15;
            wb.Title=['Calculating for ', num2str(i), '/', num2str(length(k)), ' Dt value'];
            ki=k(i); eki=ek(i);
            kn=normrnd(ki,eki/2,ttrials,1);
            Dt=10.^kn;
            for j=1:ttrials %Monte Carlo method
                if wb.CancelRequested
                    isbreak=1;
                    break
                end
                wb.Message='Monte Carlo estimating time uncertainty';
                wb.Value=j/ttrials;
                %relationship between T and t
                switch coolingpath
                    case 'isothermal'
                        D =@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-0*a)); %consider curve-fitting, temperature, and experimental errors
                    case 'linear cooling'
                        D =@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-dT_dt*a));
                    case 'exponential cooling'
                        D =@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)*exp(-dT_dt*a)));
                    case 'parabolic cooling'
                        D =@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-dT_dt*(a^2)));
                end
                tn(j) = fzero(@(a)integral(D, 0, a)-Dt(j), 1);
            end
            MC_data(:,i*2)=tn;
            tn(isnan(tn))=[];
            logtn=log(tn);
            sd=std(logtn)*2; %2 sigma absoute error
            u=mean(logtn);
            t(i)=exp(u);
            et(i,1)=(exp(u-sd)/t(i)-1)*100; %2 sigma percent error
            et(i,2)=(exp(u+sd)/t(i)-1)*100;
        end
    end
    %------------------------------------------------
    close(wb);
    if isbreak
        return
    end
end
% ==============================================================
%check whether ln(t) is normal distribution
normaltest=cell(1,length(k));
for i=1:length(k)
    tp=MC_data(:,i*2);
    tp(isnan(tp))=[];
    logtp=log(tp);
    sd=std(logtp)*2; %2 sigma absoute error
    u=mean(logtp);
    if isempty(logtp)
        normaltest(i)={'not normal distribution'};
        continue
    end
    try
        h=lillietest(logtp, 0.05);
    catch
        h=1;
    end
    if h % h = 1 indicates rejection of the null hypothesis at the 5% significance level
        normaltest(i)={'not normal distribution'};
    else % h = 0 indicates that lillietest does not reject the null hypothesis at the 5% significance level.
        normaltest(i)={'normal distribution'};
    end
    %plot figure----------------------------------
    if length(k)==1 && dT_dt~=0
        figure('Name', 'ln[t] distribution')
        histogram(logtp);
        xlabel('ln[t]')
        ylabel('Number')
        if h
            text(0.98,0.9,['ln[t] = ', num2str(u,'%.2f'), ' ± ', num2str(sd,'%.2f'), newline, ...
                '(2SD, N = ', num2str(length(logtp)), ')', newline, ...
                'not normal distribution'],...
                'Units','normalized','HorizontalAlignment',"right")
        else
            text(0.98,0.9,['ln[t] = ', num2str(u,'%.2f'), ' ± ', num2str(sd,'%.2f'), newline, ...
                '(2SD, N = ', num2str(length(logtp)), ')', newline, ...
                'normal distribution'],...
                'Units','normalized','HorizontalAlignment',"right")
        end
    end
end

%% ------------------------------------------------
name=cell(1,length(k)*2);
tstr=cell(1,length(k)*2);
for i=1:length(k)
    name{i*2-1}='Temperature (Celcius)';
    name{i*2}=['log10(Dt)=',num2str(k(i))];
    tstr{i*2-1}='';
    tstr{i*2}=['t(s), ',normaltest{i}];
end
MC_output=num2cell(MC_data);
%write into file----------------------------------
filename=['Monte Carlo result ',datestr(now,'yyyy-mm-dd HH-MM-SS.FFF AM')];
[FileName,PathName,~] = uiputfile({'*.txt';'*.csv';'*.xls';'*.xlsx'},'Save result',filename);
if ~FileName
    return
end
[~,~,extension]=fileparts(FileName);
switch extension
    case '.txt'
        writecell([name;tstr;MC_output],fullfile(PathName,FileName),"Delimiter",'tab');
    case '.csv'
        writecell([name;tstr;MC_output],fullfile(PathName,FileName));
    case {'.xlsx','.xls'}
        writecell([name;tstr;MC_output],fullfile(PathName,FileName));
        %writecell([name;tstr;MC_output],fullfile(PathName,FileName),'AutoFitWidth', 0);
        %xlswrite(fullfile(PathName,FileName),[name;tstr;MC_output]);
    otherwise
        return
end
