function [k, ek, R2] = Diffsusion_solution_AH(x, C, C1, C2, type, weights)
%% curve fit
% [x,I]=sort(x);
% C=C(I);
[fitresult, gof] = Fit_Diffusion_AH(x, C, C1, C2, type, weights);
R2=gof.rsquare;
Coef_boundary=confint(fitresult);
u=mean(Coef_boundary,1);
se=Coef_boundary(2,:)-u;
k=u(1); x0=u(2);
ek=se(1);

%% Plot fit with data.
figure;
x_p=linspace(min(x),max(x),100);
[~, C_pred] = predint(fitresult, x_p);
if find(weights~=1)
    errorbar(x, C, 1./sqrt(weights), 'LineStyle', 'none');
else
    plot(x, C, 'ko')
end
hold on
plot(x_p, C_pred, 'r-')
h=gca;
ym=h.YLim;
plot([x0,x0], ym, '--')
hold off
legend('raw data', ['curve fit, R2 = ', num2str(R2)], ['x0 = ', num2str(x0)]);
xlabel( 'X (m)' );
ylabel( 'Concentration' );
title('Curve fitting')