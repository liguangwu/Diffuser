function [k_data, ek_data, paras] = Diffsusion_solution_IL(x, C, C1, C2, C3, type, flat, weights)
%% curve fit
% [x,I]=sort(x);
% C=C(I);
if flat
    C0=C2;
else
    if strcmp(type, 'I') || strcmp(type, 'K')
        C0=max(C); %estimate initial C0, high value in the region
    else
        C0=min(C); %estimate initial C0, low value the region
    end
end
ntrials=100;
tolerance=0.001; %if >tol, accept
R2=zeros(ntrials,1);
C_data=ones(ntrials,1)*C0;
h_data=zeros(ntrials,1);
x0_data=zeros(ntrials,1);
k_data=zeros(ntrials,1);
ek_data=zeros(ntrials,1);

%estimate X0
minX=2*min(x)-max(x); maxX=max(x);
for i=1:ntrials
    x0=(maxX+minX)/2;
    x0_data(i)=x0;
    x1=x-x0;
    [fitresult, gof] = Fit_Diffusion_IL(x1, C, C1, C0, C3, type, weights);
    Coef_boundary=confint(fitresult);
    u=mean(Coef_boundary,1);
    se=Coef_boundary(2,:)-u;
    h_data(i)=u(1); k_data(i)=u(2);
    ek_data(i)=se(2);
    R2(i)=gof.rsquare;

    if i==1 %find the best fit h range using bisection method
        minX=x0;
    else
        if R2(i)>=max(R2(1:i-1)) || (sum(R2(1:i)<0)==i) %|| g1(i)<0.8
            minX=x0;
        else
            maxX=x0;
        end
    end
    
    %tolerance
    if i>1
        difference=abs(x0_data(i)-x0_data(i-1))/x0_data(i);
        if difference<tolerance
            break
        end
    end
end
R2(i+1:end)=[];
C_data(i+1:end)=[];
h_data(i+1:end)=[];
x0_data(i+1:end)=[];
k_data(i+1:end)=[];
ek_data(i+1:end)=[];
paras=[C_data, x0_data, h_data, R2];

% plot figure-------------------------------------
figure
plot(R2, '-b.')
yyaxis right
plot(x0_data, '-r.')
legend({'R2','x0'})
title('Estimate x0')
xlabel('Trials')

if flat
    %last caculation result----------------------------
    figure
    x_p=linspace(min(x1),max(x1),100);
    [~, C_pred] = predint(fitresult, x_p);
    if find(weights~=1)
        errorbar(x, C, 1./sqrt(weights), 'LineStyle', 'none');
    else
        plot(x, C, 'ko')
    end
    hold on
    plot(x_p+x0, C_pred, 'r-')
    h=gca;
    ym=h.YLim;
    plot([x0,x0], ym, '--')
    hold off
    legend('raw data', ['curve fit, R2 = ', num2str(gof.rsquare)], 'x0 position');
    xlabel( 'X (m)' );
    ylabel( 'Concentration' );
    title({['C0 = ', num2str(C0)]; ...
        ['h = ', num2str(h_data(end))]; ...
        ['x0 = ', num2str(x0)]})
    return
end
%% before estimate t
answer = questdlg('There is no flat peak or trough (C0) defined in the diffusion profile, please enter the range of C0', ...
	'Question', ...
	'Continue','Cancel','Cancel');
if strcmp(answer, 'Cancel')
    return
end
% opts.Resize='on';
opts.Interpreter='tex';
answer=inputdlg({'Min value of C0','Max value of C0', 'Step length (\DeltaC)'},'Input the range of C0', [1,50], {'','',''}, opts);
rangedata=cell2num(answer);
if isempty(answer)
    return
end
if ~isnan(rangedata(1)) && ~isnan(rangedata(2)) && ~isnan(rangedata(3))
    if strcmp(type, 'I') || strcmp(type, 'K') %high value in the region
        if rangedata(1)<C0
            errordlg(['The min value should be greater than ', num2str(C0)])
            return
        end
    else %low value in the region
        if rangedata(1)<0
            errordlg('The min value should be greater than 0')
            return
        end
        if rangedata(2)>C0
            errordlg(['The max value should be less than ', num2str(C0)])
            return
        end
    end
    if rangedata(3)<0
            errordlg('The step length should be greater than 0')
            return
    end
else
    errordlg('Please input a proper concentration range and step length')
    return
end

if strcmp(type, 'I') || strcmp(type, 'K') %high value in the region
    C_data=rangedata(1):rangedata(3):rangedata(2);
else %low value in the region
    C_data=rangedata(2):-rangedata(3):rangedata(1);
end
C_data=C_data';
%% estimate t
ntrials=length(C_data);
R2=zeros(ntrials,1);
h_data=zeros(ntrials,1);
x0_data=ones(ntrials,1)*x0; %fixed at this step
k_data=zeros(ntrials,1);
ek_data=zeros(ntrials,1);

for i=1:ntrials
    C0=C_data(i);
    [fitresult, gof] = Fit_Diffusion_IL(x1, C, C1, C0, C3, type, weights);
    Coef_boundary=confint(fitresult);
    u=mean(Coef_boundary,1);
    se=Coef_boundary(2,:)-u;
    h_data(i)=u(1); k_data(i)=u(2);
    ek_data(i)=se(2);
    R2(i)=gof.rsquare;
end

% plot figure-------------------------------------
figure
plot(C_data, R2, '-b.');
yyaxis right
plot(C_data, k_data, '-m.');
legend({'R2','k'})
yyaxis left
title(['Estiamte k when x0 is fixed at ', num2str(x0)])
xlabel('C0')

% last caculation result----------------------------
figure
x_p=linspace(min(x1),max(x1),100);
[~, C_pred] = predint(fitresult,x_p);
if find(weights~=1)
    errorbar(x, C, 1./sqrt(weights), 'LineStyle', 'none');
else
    plot(x, C, 'ko')
end
hold on
plot(x_p+x0, C_pred, 'r-')
h=gca;
ym=h.YLim;
plot([x0,x0], ym, '--')
hold off
legend('raw data', ['curve fit, R2 = ', num2str(gof.rsquare)], 'x0 position');
xlabel( 'X (m)' );
ylabel( 'Concentration' );
title({['last trial, C0 = ', num2str(C0)]; ...
    ['h = ', num2str(h_data(end))]; ...
    ['x0 = ', num2str(x0)]})
paras=[C_data, x0_data, h_data, R2];