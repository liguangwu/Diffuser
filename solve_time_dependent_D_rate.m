function [timescale, et, MC_t, rate, erate, MC_rate, error_contribution, isbreak]= solve_time_dependent_D_rate(k, ek, initialT, eT, ttrials, DiffCoef, coolingpath, withD0Eerror, app)
%calculate rate automatically for non-isothermal condition
%ek, and eT are 2sigma absolute errors, et and erate are 2 sigma percent errors [-2s%, +2s%]
%k=log10[Dt]
isbreak=0; MC_t=[]; MC_rate=[]; error_contribution=[];
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
Tn1=normrnd(initialT,0,ttrials,1);
Tn2=normrnd(initialT,eT/2,ttrials,1);
R1 = mvnrnd([lnD0,E],cov1,ttrials);
D0n1=exp(R1(:,1)); En1=R1(:,2); 
R2 = mvnrnd([lnD0,E],cov2,ttrials);
D0n2=exp(R2(:,1)); En2=R2(:,2);
timescale=NaN(length(k),1);
et=NaN(length(k),2);
MC_t=NaN(ttrials, length(k)*2);
rate=NaN(length(k),1);
erate=NaN(length(k),2);
MC_rate=NaN(ttrials, length(k)*2);
% ==============================================================
tn=NaN(ttrials,1); tn1=NaN(ttrials,1); tn2=NaN(ttrials,1); tn3=NaN(ttrials,1);
raten=NaN(ttrials,1); raten1=NaN(ttrials,1); raten2=NaN(ttrials,1); raten3=NaN(ttrials,1);
wb=uiprogressdlg(app.DiffuserGUI,'ShowPercentage', 'on', 'Cancelable', 'on');

if length(k)==1
    i=1;
    MC_t(:,i*2-1)=Tn2-273.15;
    MC_rate(:,i*2-1)=Tn2-273.15;
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
            case 'linear cooling'
                raten1(j)=fminsearch(@(v) abs(log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-v*a)), 0, (Tn1(j)-273)/v))-log(Dt(j))), 1); %tmax=T0/v
                raten2(j)=fminsearch(@(v) abs(log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-v*a)), 0, (Tn2(j)-273)/v))-log(Dt(j))), 1);
                raten3(j)=fminsearch(@(v) abs(log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-v*a)), 0, (Tn2(j)-273)/v))-log(Dt(j))), 1);
                %tmp=integral(fun,0,tmax), so tmp-integral(fun,0,t) must>=0.
                %To minimize tmp-integral(fun,0,t) and find the minumum t:
                % tn1(j) = fminbnd(@(t) -log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-raten1(j)*a)), 0, t))+log(Dt(j)), 0, (Tn1(j)-273)/raten1(j));
                % tn2(j) = fminbnd(@(t) -log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-raten2(j)*a)), 0, t))+log(Dt(j)), 0, (Tn2(j)-273)/raten2(j));
                % tn3(j) = fminbnd(@(t) -log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-raten3(j)*a)), 0, t))+log(Dt(j)), 0, (Tn2(j)-273)/raten3(j));
                tmp1=log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-raten1(j)*a)), 0, (Tn1(j)-273)/raten1(j)));
                tn1(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-raten1(j)*a)), 0, t))+tmp1), 0, (Tn1(j)-273)/raten1(j));
                tmp2=log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-raten2(j)*a)), 0, (Tn2(j)-273)/raten2(j)));
                tn2(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-raten2(j)*a)), 0, t))+tmp2), 0, (Tn2(j)-273)/raten2(j));
                tmp3=log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-raten3(j)*a)), 0, (Tn2(j)-273)/raten3(j)));
                tn3(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-raten3(j)*a)), 0, t))+tmp3), 0, (Tn2(j)-273)/raten3(j));              
            case 'exponential cooling'
                raten1(j)=fminsearch(@(v) abs(log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)*exp(-v*a))), 0, log(Tn1(j)/273)/v))-log(Dt(j))), 1); %tmax=ln(T0)/v
                raten2(j)=fminsearch(@(v) abs(log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)*exp(-v*a))), 0, log(Tn2(j)/273)/v))-log(Dt(j))), 1);
                raten3(j)=fminsearch(@(v) abs(log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)*exp(-v*a))), 0, log(Tn2(j)/273)/v))-log(Dt(j))), 1);
                % tn1(j)=fminbnd(@(t,a) -log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)*exp(-raten1(j)*a))), 0, t))+log(Dt(j)), 0, log(Tn1(j)/273)/raten1(j));
                % tn2(j)=fminbnd(@(t,a) -log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)*exp(-raten2(j)*a))), 0, t))+log(Dt(j)), 0, log(Tn2(j)/273)/raten2(j));
                % tn3(j)=fminbnd(@(t,a) -log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)*exp(-raten3(j)*a))), 0, t))+log(Dt(j)), 0, log(Tn2(j)/273)/raten3(j));
                tmp1=log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)*exp(-raten1(j)*a))), 0, log(Tn1(j)/273)/raten1(j)));
                tn1(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)*exp(-raten1(j)*a))), 0, t))+tmp1), 0, log(Tn1(j)/273)/raten1(j));
                tmp2=log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)*exp(-raten2(j)*a))), 0, log(Tn2(j)/273)/raten2(j)));
                tn2(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)*exp(-raten2(j)*a))), 0, t))+tmp2), 0, log(Tn2(j)/273)/raten2(j));
                tmp3=log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)*exp(-raten3(j)*a))), 0, log(Tn2(j)/273)/raten3(j)));
                tn3(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)*exp(-raten3(j)*a))), 0, t))+tmp3), 0, log(Tn2(j)/273)/raten3(j)); 
            case 'parabolic cooling'
                raten1(j)=fminsearch(@(v) abs(log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-v*(a.^2))), 0, sqrt((Tn1(j)-273)/v)))-log(Dt(j))), 1); %tmax=sqrt(T0/v)
                raten2(j)=fminsearch(@(v) abs(log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-v*(a.^2))), 0, sqrt((Tn2(j)-273)/v)))-log(Dt(j))), 1);
                raten3(j)=fminsearch(@(v) abs(log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-v*(a.^2))), 0, sqrt((Tn2(j)-273)/v)))-log(Dt(j))), 1);
                % tn1(j)=fminbnd(@(t,a) -log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-raten1(j)*(a.^2))), 0, t))+log(Dt(j)), 0, sqrt((Tn1(j)-273)/raten1(j)));
                % tn2(j)=fminbnd(@(t,a) -log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-raten2(j)*(a.^2))), 0, t))+log(Dt(j)), 0, sqrt((Tn2(j)-273)/raten2(j)));
                % tn3(j)=fminbnd(@(t,a) -log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-raten3(j)*(a.^2))), 0, t))+log(Dt(j)), 0, sqrt((Tn2(j)-273)/raten3(j)));
                tmp1=log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-raten1(j)*(a.^2))), 0, sqrt((Tn1(j)-273)/raten1(j))));
                tn1(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn1(j)-raten1(j)*(a.^2))), 0, t))+tmp1), 0, sqrt((Tn1(j)-273)/raten1(j)));
                tmp2=log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-raten2(j)*(a.^2))), 0, sqrt((Tn2(j)-273)/raten2(j))));
                tn2(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n1(j)*exp(-En1(j)*1000/8.314./(Tn2(j)-raten2(j)*(a.^2))), 0, t))+tmp2), 0, sqrt((Tn2(j)-273)/raten2(j)));
                tmp3=log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-raten3(j)*(a.^2))), 0, sqrt((Tn2(j)-273)/raten3(j))));
                tn3(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n2(j)*exp(-En2(j)*1000/8.314./(Tn2(j)-raten3(j)*(a.^2))), 0, t))+tmp3), 0, sqrt((Tn2(j)-273)/raten3(j))); 
        end
    end
    %error contribution------------------------------------------------
    if withD0Eerror
        %timescale
        MC_t(:,i*2)=tn3;
        tn1(isnan(tn1))=[]; tn2(isnan(tn2))=[]; tn3(isnan(tn3))=[];
        logtn1=log(tn1); logtn2=log(tn2); logtn3=log(tn3);
        sd1=std(logtn1)*2; sd2=std(logtn2)*2; sd3=std(logtn3)*2; %2 sigma absoute error
        u1=mean(logtn1); u2=mean(logtn2); u3=mean(logtn3);
        timescale(i)=exp(u3);
        et(i,1)=(exp(u3-sd3)/timescale(i)-1)*100; %2 sigma percent error
        et(i,2)=(exp(u3+sd3)/timescale(i)-1)*100;
        %rate
        MC_rate(:,i*2)=raten3;
        raten1(isnan(raten1))=[]; raten2(isnan(raten2))=[]; raten3(isnan(raten3))=[];
        lograten1=log(raten1); lograten2=log(raten2); lograten3=log(raten3);
        sd1_r=std(lograten1)*2; sd2_r=std(lograten2)*2; sd3_r=std(lograten3)*2; %2 sigma absoute error
        u1_r=mean(lograten1); u2_r=mean(lograten2); u3_r=mean(lograten3);
        rate(i)=exp(u3_r);
        erate(i,1)=(exp(u3_r-sd3_r)/rate(i)-1)*100; %2 sigma percent error
        erate(i,2)=(exp(u3_r+sd3_r)/rate(i)-1)*100;
    else
        %timescale
        MC_t(:,i*2)=tn2;
        tn1(isnan(tn1))=[]; tn2(isnan(tn2))=[];
        logtn1=log(tn1); logtn2=log(tn2);
        sd1=std(logtn1)*2; sd2=std(logtn2)*2; %2 sigma absoute error
        u1=mean(logtn1); u2=mean(logtn2);
        timescale(i)=exp(u2);
        et(i,1)=(exp(u2-sd2)/timescale(i)-1)*100; %2 sigma percent error
        et(i,2)=(exp(u2+sd2)/timescale(i)-1)*100;
        u3=u2; sd3=sd2;
        %rate
        MC_rate(:,i*2)=raten2;
        raten1(isnan(raten1))=[]; raten2(isnan(raten2))=[];
        lograten1=log(raten1); lograten2=log(raten2);
        sd1_r=std(lograten1)*2; sd2_r=std(lograten2)*2; %2 sigma absoute error
        u1_r=mean(lograten1); u2_r=mean(lograten2);
        rate(i)=exp(u2_r);
        erate(i,1)=(exp(u2_r-sd2_r)/rate(i)-1)*100; %2 sigma percent error
        erate(i,2)=(exp(u2_r+sd2_r)/rate(i)-1)*100;
        u3_r=u2_r; sd3_r=sd2_r;
    end
    sd1=sd1+u3-u1; sd2=sd2+u3-u2; %keep the average same
    if sd1<=sd2 && sd2<=sd3 %otherwise, the trials are not enough
        e1=sd1/sd3*100; e2=(sd2/sd3-sd1/sd3)*100; e3=(1-sd2/sd3)*100; %contribution to ln(t) error from curve fitting, temperature, and experimental errors
        error_contribution(1,:)=[e1,e2,e3];
    end
    sd1_r=sd1_r+u3_r-u1_r; sd2_r=sd2_r+u3_r-u2_r; %keep the average same
    if sd1_r<=sd2_r && sd2_r<=sd3_r %otherwise, the trials are not enough
        e1=sd1_r/sd3_r*100; e2=(sd2_r/sd3_r-sd1_r/sd3_r)*100; e3=(1-sd2_r/sd3_r)*100; %contribution to ln(t) error from curve fitting, temperature, and experimental errors
        error_contribution(2,:)=[e1,e2,e3];
    end
end
%------------------------------------------------
%calculate t when no flat peak/trough in profiles I-L
%consider curve-fitting, temperature, and experimental errors if selected,
%but won't calculate error contribution
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
        MC_t(:,i*2-1)=Tn-273.15;
        MC_rate(:,i*2-1)=Tn-273.15;
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
                case 'linear cooling'
                    raten(j)=fminsearch(@(v) abs(log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-v*a)), 0, (Tn(j)-273)/v))-log(Dt(j))), 1); %tmax=T0/v
                    tmp=log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-raten(j)*a)), 0, (Tn(j)-273)/raten(j)));
                    tn(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-raten(j)*a)), 0, t))+tmp), 0, (Tn(j)-273)/raten(j));
                case 'exponential cooling'
                    raten(j)=fminsearch(@(v) abs(log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)*exp(-v*a))), 0, log(Tn(j)/273)/v))-log(Dt(j))), 1); %tmax=ln(T0)/v
                    tmp=log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)*exp(-raten(j)*a))), 0, log(Tn(j)/273)/raten(j)));
                    tn(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)*exp(-raten(j)*a))), 0, t))+tmp), 0, log(Tn(j)/273)/raten(j));
                case 'parabolic cooling'
                    raten(j)=fminsearch(@(v) abs(log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-v*(a.^2))), 0, sqrt(Tn1(j)/v)))-log(Dt(j))), 1); %tmax=sqrt(T0/v)
                    tmp=log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-raten(j)*(a.^2))), 0, sqrt((Tn(j)-273)/raten(j))));
                    tn(j) = fminbnd(@(t) abs(-log(integral(@(a) D0n(j)*exp(-En(j)*1000/8.314./(Tn(j)-raten(j)*(a.^2))), 0, t))+tmp), 0, sqrt((Tn(j)-273)/raten(j)));
            end
        end
        MC_t(:,i*2)=tn; MC_rate(:,i*2)=raten;
        tn(isnan(tn))=[]; raten(isnan(raten))=[];
        logtn=log(tn); lograten=log(raten);
        sd=std(logtn)*2; sd_r=std(lograten)*2;%2 sigma absoute error
        u=mean(logtn); u_r=mean(lograten);
        timescale(i)=exp(u);
        et(i,1)=(exp(u-sd)/timescale(i)-1)*100; %2 sigma percent error
        et(i,2)=(exp(u+sd)/timescale(i)-1)*100;
        rate(i)=exp(u_r);
        erate(i,1)=(exp(u_r-sd_r)/rate(i)-1)*100; %2 sigma percent error
        erate(i,2)=(exp(u_r+sd_r)/rate(i)-1)*100;
    end
end
%------------------------------------------------
close(wb);
if isbreak
    return
end

% ==============================================================
%check whether ln(t) is normal distribution
normaltest=cell(1,length(k));
normaltest_r=cell(1,length(k));
for i=1:length(k)
    tp=MC_t(:,i*2); tp_r=MC_rate(:,i*2);
    tp(isnan(tp))=[]; tp_r(isnan(tp_r))=[];
    logtp=log(tp); logtp_r=log(tp_r);
    sd=std(logtp)*2; sd_r=std(logtp_r)*2; %2 sigma absoute error
    u=mean(logtp); u_r=mean(logtp_r);
    if isempty(logtp)
        normaltest(i)={'not normal distribution'};
        continue
    end
    if isempty(logtp_r)
        normaltest_r(i)={'not normal distribution'};
        continue
    end
    try
        h=lillietest(logtp, 0.05);
        h_r=lillietest(logtp_r, 0.05);
    catch
        h=1;
        h_r=1;
    end
    if h % h = 1 indicates rejection of the null hypothesis at the 5% significance level
        normaltest(i)={'not normal distribution'};
    else % h = 0 indicates that lillietest does not reject the null hypothesis at the 5% significance level.
        normaltest(i)={'normal distribution'};
    end
    if h_r % h = 1 indicates rejection of the null hypothesis at the 5% significance level
        normaltest_r(i)={'not normal distribution'};
    else % h = 0 indicates that lillietest does not reject the null hypothesis at the 5% significance level.
        normaltest_r(i)={'normal distribution'};
    end
    %plot figure----------------------------------
    if length(k)==1
        %check timescale
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
        %check cooling rate
        figure('Name', 'ln[rate] distribution')
        histogram(logtp_r);
        xlabel('ln[rate]')
        ylabel('Number')
        if h
            text(0.98,0.9,['ln[rate] = ', num2str(u_r,'%.2f'), ' ± ', num2str(sd_r,'%.2f'), newline, ...
                '(2SD, N = ', num2str(length(logtp_r)), ')', newline, ...
                'not normal distribution'],...
                'Units','normalized','HorizontalAlignment',"right")
        else
            text(0.98,0.9,['ln[rate] = ', num2str(u_r,'%.2f'), ' ± ', num2str(sd_r,'%.2f'), newline, ...
                '(2SD, N = ', num2str(length(logtp_r)), ')', newline, ...
                'normal distribution'],...
                'Units','normalized','HorizontalAlignment',"right")
        end
    end
end

%% ------------------------------------------------
name=cell(1,length(k)*3);
tstr=cell(1,length(k)*3);
MC_output=[];
for i=1:length(k)
    name{i*3-2}=['log10(Dt)=',num2str(k(i))];
    name{i*3-1}='';
    name{i*3}='';
    tstr{i*3-2}='Temperature (Celcius)';
    tstr{i*3-1}=['t(s), ',normaltest{i}];
    switch coolingpath
        case 'linear cooling'
            tstr{i*3}=['Cooling rate (Celcius/s), ',normaltest_r{i}];
        case 'exponential cooling'
            tstr{i*3}=['Cooling rate (1/s), ',normaltest_r{i}];
        case 'parabolic cooling'
            tstr{i*3}=['Cooling rate (Celcius/s^2), ',normaltest_r{i}];
    end
    MC_output=cat(2, MC_output, num2cell([MC_t(:,(i*2-1):(i*2)), MC_rate(:,i*2)]));
end

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
