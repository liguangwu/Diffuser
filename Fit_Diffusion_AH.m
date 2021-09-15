function [fitresult, gof] = Fit_Diffusion_AH(x, C, C1, C2, type, weights)

[xData, yData] = prepareCurveData( x, C );
% Set up fittype and options.
%half profile
%     C =(C2-C1)*erf(x/k) + C1; %x>0, Crim<Ccore
%     C =(C2-C1)*erf(-x/k) + C1; %x<0, Crim<Ccore
%     C =(C2-C1)*erfc(x/k) + C1; %x>0, Crim>Ccore
%     C =(C2-C1)*erfc(-x/k) + C1; %x<0, Crim>Ccore
%full profile
%     C =(C2-C1)*(erf((h+x)/k) + erf((h-x)/k)-1) + C1; %Crim<Ccore
%     C =(C2-C1)*(erfc((h+x)/k) + erfc((h-x)/k)) + C1; %Crim>Ccore
h=(max(x)-min(x))/2;
x0=(max(x)+min(x))/2; %initial x0 guess
k0=1e-7; %initial k guess

switch type
    case 'A'
        ft = fittype( [num2str(abs(C2-C1)),'*erf((x-x0)/k)+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    case 'B'
        ft = fittype( [num2str(abs(C2-C1)),'*erf((-x+x0)/k)+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    case 'C'
        ft = fittype( [num2str(abs(C2-C1)),'*erfc((x-x0)/k)+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    case 'D'
        ft = fittype( [num2str(abs(C2-C1)),'*erfc((-x+x0)/k)+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    case 'E'
        ft = fittype( [num2str(abs(C2-C1)),'*(erf((', num2str(h), '+x-x0)/k)+erf((', num2str(h), '-x+x0)/k)-1)+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    case 'F'
        ft = fittype( [num2str(abs(C2-C1)),'*(erfc((', num2str(h), '+x-x0)/k)+erfc((', num2str(h), '-x+x0)/k))+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    case 'G'
        ft = fittype( [num2str(abs((C2-C1))/2),'*(1+erf((x-x0)/k))+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    case 'H'
        ft = fittype( [num2str(abs((C1-C2))/2),'*erfc((x-x0)/k)+', num2str(min(C1,C2))], 'independent', 'x', 'dependent', 'y' );
    otherwise
        return
end

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [k0 x0];
opts.Weights = weights;
opts.Lower=[0,-Inf];
opts.Upper=[Inf,Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure;
% h = plot( fitresult, xData, yData );
% legend( h, 'raw data', 'curve fit');
% xlabel( 'x' );
% ylabel( 'C' );
