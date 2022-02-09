function [log10Dt, elog10Dt, paras] = solution_AH(x, C, C1, C2, profiletype, weights, factor, app, model_x, props_matrix, pmsum)
%% curve fit
% [x,I]=sort(x);
% C=C(I);
% weights=weights(I);
if isempty(model_x)
    [p, ci, C_pred, xp, yp] = Fit_Diffusion_AH(x, C, C1, C2, profiletype, weights, factor);
else
    wb=uiprogressdlg(app.DiffuserGUI,'Title','Deconvoluting the diffusion profile','Indeterminate','on');
    [p, ci, C_pred, yp] = deconv_AH(x, C, C1, C2, profiletype, weights, factor, model_x, props_matrix, pmsum);
    close(wb)
    xp=model_x;
    yp(xp<min(x),:)=[];
    xp(xp<min(x))=[];
    yp(xp>max(x),:)=[];
    xp(xp>max(x))=[];
end
sse=sum(weights.*((C-C_pred).^2));
sst=sum(weights.*(C-mean(C)).^2);
R2=1-sse/sst;
se=ci(:,2)-p';
log10Dt=p(4); 
elog10Dt=se(4)/tinv(0.975,length(x)-length(p)+length([C1,C2]))*2;
paras=[R2,p]; %paras=R2; C1; C2; x0; log10(Dt); h
%% Plot fit with data.
figure;
if find(weights~=1)
    h=errorbar(x, C, 1./sqrt(weights));
    h.Marker='.';
    h.MarkerSize=12;
    h.LineStyle='none';
    h.CapSize=0;
else
    plot(x, C, 'ko')
end
hold on
%initial concentration----------------------------
%p=C1; C2; x0; log10(Dt); h
C1=p(1); C2=p(2); x0=p(3); 
switch profiletype
    case 'A'
        plot([x0,x0,max(x)],[C1,C2,C2],'k--');
    case 'B'
        plot([x0,x0,min(x)],[C1,C2,C2],'k--');
    case 'C'
        plot([x0,x0,max(x)],[C2,C1,C1],'k--');
    case 'D'
        plot([x0,x0,min(x)],[C2,C1,C1],'k--');
    case 'E'
        h=p(5);
        plot([x0-h,x0-h,x0+h,x0+h],[C1,C2,C2,C1],'k--');
    case 'F'
        h=p(5);
        plot([x0-h,x0-h,x0+h,x0+h],[C2,C1,C1,C2],'k--');
    case 'G'
        plot([min(x),x0,x0,max(x)],[C1,C1,C2,C2],'k--');
    case 'H'
        plot([min(x),x0,x0,max(x)],[C2,C2,C1,C1],'k--');
end
%fitted curve------------------------------------
hp=plot(xp, yp, 'm');
hold off

if ~isempty(model_x)
    hp(1).Color='b'; hp(1).LineStyle='-';
    legend('Data', 'Initial condition', 'Convoluted', ['Deconvoluted:', newline, 'log_{10}[Dt] = ', num2str(log10Dt), newline, 'R^{2} = ', num2str(R2)], ...
        'Location','best');
    if 2*sqrt(10^log10Dt)>mean(diff(x))
        title('Curve fitting with deconvolution')
    else
        title('Curve fitting with deconvolution: no resolvable profile')
    end
else
    legend('Data', 'Initial condition', ['Curve fit:', newline, 'log_{10}[Dt] = ', num2str(log10Dt), newline, 'R^{2} = ', num2str(R2)], ...
        'Location','best');
    title('Curve fitting')
end
xlabel( 'X (m)' );
ylabel( 'Concentration' );
