function [fitresult, gof] = Fit_Diffusion_IL(x, C, C1, C0, C3, type, weights)

[xData, yData] = prepareCurveData( x, C );
h0=(max(x)-min(x))/2; %initial h guess
% Set up fittype and options.
switch type
    case 'I'
        ft = fittype( [num2str((C0-C3)/2), '*erf((h-x)/k)+', num2str((C0-C1)/2), '*erf((h+x)/k)+', ...
            num2str((C1+C3)/2)], 'independent', 'x', 'dependent', 'y' );
        k0=h0/erfinv( (mean(C)-(C1+C3)/2) / (abs(C0-C3)/2+abs(C0-C1)/2) ); %initial k0 guess, erfinv(A), |A|<1
    case 'J'
        ft = fittype( [num2str(abs(C0-C3)/2), '*erfc((h-x)/k)+', num2str(abs(C0-C1)/2), '*erfc((h+x)/k)+', ...
            num2str(C0)], 'independent', 'x', 'dependent', 'y' );
        k0=h0/erfcinv( (mean(C)-C0) / (abs(C0-C3)/2+abs(C0-C1)/2) ); %initial k0 guess
    case 'K'
        ft = fittype( [num2str((C0-C3)/2), '*erf((h-x)/k)+', num2str((C0-C1)/2), '*erf((h+x)/k)+', ...
            num2str((C1+C3)/2)], 'independent', 'x', 'dependent', 'y' );
        k0=h0/erfinv( (mean(C)-(C1+C3)/2) / (abs(C0-C3)/2+abs(C0-C1)/2) ); %initial k0 guess
    case 'L'
        ft = fittype( [num2str(abs(C0-C3)/2), '*erfc((h-x)/k)+', num2str(abs(C0-C1)/2), '*erfc((h+x)/k)+', ...
            num2str(C0)], 'independent', 'x', 'dependent', 'y' );
        k0=h0/erfcinv( (mean(C)-C0) / (abs(C0-C3)/2+abs(C0-C1)/2) ); %initial k0 guess
    otherwise
        return
end
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [h0, k0];
opts.Weights = weights;
opts.Lower=[0,0];
opts.Upper=[max(x)-min(x),Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure;
% h = plot( fitresult, xData, yData );
% legend( h, 'raw data', 'curve fit');
% xlabel( 'x' );
% ylabel( 'C' );


