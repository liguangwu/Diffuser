function [log10Dt, elog10Dt, paras] = solution_IL(x, C, C1, C2, C3, profiletype, weights, factor, app, model_x, props_matrix, pmsum)
%% curve fit
% [x,I]=sort(x);
% C=C(I);
% weights=weights(I);
log10Dt=[]; elog10Dt=[]; paras=[];
flat=app.AflatpeaktroughcanbeseeninprofileILCheckBox.Value;
if ~flat
    rangedata=[app.minEditField.Value, app.maxEditField.Value, app.steplengthEditField.Value];
    C_data=rangedata(1):rangedata(3):rangedata(2);
    C_data=C_data';
    if isempty(C_data)
        uialert(app.DiffuserGUI, 'Please input an appropriate concentration range for the initial flat peak/trough','');
        return
    end
    if contains(profiletype, {'I','K'}) %high value in the region
        C2=C_data;
    else %low value in the region
        C1=C_data;
    end
end
R2=[]; p_data=[];
if flat
    if isempty(model_x)
        [p, ci, C_pred, xp, yp] = Fit_Diffusion_IL(x, C, C1, C2, C3, profiletype, weights, factor);
    else
        wb=uiprogressdlg(app.DiffuserGUI,'Title','Deconvoluting the diffusion profile','Indeterminate','on');
        [p, ci, C_pred, yp] = deconv_IL(x, C, C1, C2, C3, profiletype, weights, factor, model_x, props_matrix, pmsum);
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
    elog10Dt=se(4)/tinv(0.975,length(x)-length(p)+length([C1,C2,C3]))*2;
    paras=[R2,p]; %paras=R2; C1; C2; x0; log10(Dt); h; C3
else
    for i=1:length(C_data)
        if contains(profiletype, {'I','K'}) %high value in the region 
            if isempty(model_x)
                [p, ci, C_pred, xp, yp] = Fit_Diffusion_IL(x, C, C1, C2(i), C3, profiletype, weights, factor);
            else
                wb=uiprogressdlg(app.DiffuserGUI,'Title','Deconvoluting the diffusion profile','Indeterminate','on');
                [p, ci, C_pred, yp] = deconv_IL(x, C, C1, C2(i), profiletype, weights, factor, model_x, props_matrix, pmsum);
                close(wb)
            end
            length_p=length(p)-length([C1, C2(i), C3]);
        else %low value in the region
            if isempty(model_x)
                wb=uiprogressdlg(app.DiffuserGUI,'Title','Deconvoluting the diffusion profile','Indeterminate','on');
                [p, ci, C_pred, xp, yp] = Fit_Diffusion_IL(x, C, C1(i), C2, C3, profiletype, weights, factor, model_x, props_matrix, pmsum);
                close(wb)
            else
                wb=uiprogressdlg(app.DiffuserGUI,'Title','Deconvoluting the diffusion profile','Indeterminate','on');
                [p, ci, C_pred, yp] = deconv_IL(x, C,  C1(i), C2, C3, profiletype, weights, factor, model_x, props_matrix, pmsum);
                close(wb)
            end
            length_p=length(p)-length([C1(i), C2, C3]);
        end 
        sse=sum(weights.*((C-C_pred).^2));
        sst=sum(weights.*(C-mean(C)).^2);
        R2=cat(1,R2,1-sse/sst);
        se=ci(:,2)-p';
        log10Dt=cat(1,log10Dt,p(4));
        elog10Dt=cat(1,elog10Dt,se(4)/tinv(0.975,length(x)-length_p)*2);
        p_data=cat(1,p_data,p);
    end
    paras=[R2,p_data]; %paras=R2; C1; C2; x0; log10(Dt); h; C3
end

%% Plot fit with data.
figure;
if flat
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
    %p= C1; C2; x0; log10(Dt); h; C3
    C1=p(1); C2=p(2); x0=p(3); h=p(5);
    switch profiletype
    case 'I'
        plot([min(x),x0-h,x0-h,x0+h,x0+h,max(x)],[C1,C1,C2,C2,C1,C1],'k--');
    case 'J'
        plot([min(x),x0-h,x0-h,x0+h,x0+h,max(x)],[C2,C2,C1,C1,C2,C2],'k--');
    case 'K'
        C3=p(6);
        if C(1)<C(end)
            plot([min(x),x0-h,x0-h,x0+h,x0+h,max(x)],[C1,C1,C2,C2,C3,C3],'k--');
        else
            plot([min(x),x0-h,x0-h,x0+h,x0+h,max(x)],[C3,C3,C2,C2,C1,C1],'k--');
        end
    case 'L'
        C3=p(6);
        if C(1)<C(end)
            plot([min(x),x0-h,x0-h,x0+h,x0+h,max(x)],[C3,C3,C1,C1,C2,C2],'k--');
        else
            plot([min(x),x0-h,x0-h,x0+h,x0+h,max(x)],[C2,C2,C1,C1,C3,C3],'k--');
        end
    end
    %fitted curve------------------------------------
    hp=plot(xp, yp, 'm');
    hold off
    
    if ~isempty(model_x)
        hp(1).Color='b'; hp(1).LineStyle='-';
        legend('Data', 'Initial condition', 'Convoluted', ['Deconvoluted:', newline, 'log_{10}(Dt) = ', num2str(log10Dt), newline, 'R^{2} = ', num2str(R2)], ...
            'Location','best');
        if 2*sqrt(10^log10Dt)>mean(diff(x))
            title('Curve fitting with deconvolution')
        else
            title('Curve fitting with deconvolution: no resolvable profile')
        end
    else
        legend('Data', 'Initial condition', ['Curve fit:', newline, 'log_{10}(Dt) = ', num2str(log10Dt), newline, 'R^{2} = ', num2str(R2)], ...
            'Location','best');
        title('Curve fitting')
    end
    xlabel( 'X (m)' );
    ylabel( 'Concentration' );
else
    plot(C_data, log10Dt, '-m.', 'MarkerSize',16);
    ylabel('log_{10}(Dt)')
    yyaxis right
    plot(C_data, R2, '-b.', 'MarkerSize',16);
    ylabel('R^{2}')
    set(gca,'Ycolor','b')
    legend({'log_{10}(Dt)','R^{2}'})
    title('Estiamte Dt with no flat peak/trough in profile I-L')
    xlabel('C0')
    yyaxis left
    set(gca,'Ycolor','m')
end
