function [p_min, ci_min, C_pred_min, model_y] = deconv_AH(x, C, C1, C2, profiletype, weights, factor, model_x, props_matrix, pmsum)
% Set up fittype and options.
% C1=p(1); C2=p(2); x0=p(3); log10(Dt)=p(4); h=p(5);
p_min=[]; ci_min=[]; C_pred_min=[]; model_y=[]; rsn_min=1e10;
%initial guess
h0=(max(x)-min(x))/2;
if contains(profiletype, {'A','C'})
    x0=min(x); lx0=-inf; ux0=min(x);
elseif contains(profiletype, {'B','D'})
    x0=max(x); lx0=max(x); ux0=inf;
else
    x0=(max(x)+min(x))/2; lx0=min(x); ux0=max(x);
end
switch factor
    case 1
        guess=2;
    case 1e3
        guess=-4;
    case 1e6
        guess=-10;
    case 1e9
        guess=-16;
end
gg=-15:15;
warning('off')
opts = optimset('Display','off');
% 1. unknown C1 and C2==============================================================
if isempty(C1) && isempty(C2)
    fun = @(p) sqrt(weights).*(f_convolution1(x,model_x,props_matrix,pmsum,profiletype,p)-C);
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'E','F'})
            p0 = [min(C),max(C),x0,guess+g,h0];
            lb=[0, min(C), lx0, -50, 0];
            ub=[max(C), inf, ux0, 20, inf];
        else
            p0 = [min(C),max(C),x0,guess+g];
            lb=[0, min(C), lx0, -50];
            ub=[max(C), inf, ux0, 20];
        end
        try
            [p,resnorm,R,~,~,~,J] = lsqnonlin(fun,p0,lb,ub,opts);
            ci = nlparci(p,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=p;
                ci_min=ci;
                C_pred_min = f_profile1(x,p,profiletype);
                [~,model_y(:,1),model_y(:,2)] = f_convolution1(x,model_x,props_matrix,pmsum,profiletype,p);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution1(x,model_x,props_matrix,pmsum,profiletype,p)
        original=f_profile1(model_x,p,profiletype);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile1(x,p,profiletype)
        switch profiletype
            case 'A'
                C = (p(2)-p(1))*erf((x-p(3))/2/sqrt(10^p(4)))+p(1);
            case 'B'
                C = (p(2)-p(1))*erf((-x+p(3))/2/sqrt(10^p(4)))+p(1);
            case 'C'
                C = (p(2)-p(1))*erfc((x-p(3))/2/sqrt(10^p(4)))+p(1);
            case 'D'
                C = (p(2)-p(1))*erfc((-x+p(3))/2/sqrt(10^p(4)))+p(1);
            case 'E'
                C = (p(2)-p(1))*(erf((p(5)+x-p(3))/2/sqrt(10^p(4)))+erf((p(5)-x+p(3))/2/sqrt(10^p(4)))-1)+p(1);
            case 'F'
                C = (p(2)-p(1))*(erfc((p(5)+x-p(3))/2/sqrt(10^p(4)))+erfc((p(5)-x+p(3))/2/sqrt(10^p(4))))+p(1);
            case 'G'
                C = (p(2)-p(1))/2*(1+erf((x-p(3))/2/sqrt(10^p(4))))+p(1);
            case 'H'
                C = (p(2)-p(1))/2*erfc((x-p(3))/2/sqrt(10^p(4)))+p(1);
        end
    end
% 2. only know C1 ==============================================================
% C2=q(1); x0=q(2); log10(Dt)=q(3); h=q(4);
if ~isempty(C1) && isempty(C2)
    fun = @(q) sqrt(weights).*(f_convolution2(x,model_x,props_matrix,pmsum,profiletype,q,C1)-C);
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'E','F'})
            q0 = [max(C),x0,guess+g,h0];
            lb=[min(C), min(x), -50, 0];
            ub=[inf, max(x), 20, inf];
        else
            q0 = [max(C),x0,guess+g];
            lb=[min(C), min(x), -50];
            ub=[inf, max(x), 20];
        end
        try
            [q,resnorm,R,~,~,~,J] = lsqnonlin(fun,q0,lb,ub,opts);
            ci = nlparci(q,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=[C1,q];
                ci_min=[C1,C1;ci];
                C_pred_min = f_profile2(x,q,profiletype,C1);
                [~,model_y(:,1),model_y(:,2)] = f_convolution2(x,model_x,props_matrix,pmsum,profiletype,q,C1);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution2(x,model_x,props_matrix,pmsum,profiletype,p,C1)
        original=f_profile2(model_x,p,profiletype,C1);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile2(x,q,profiletype,C1)
        switch profiletype
            case 'A'
                C = (q(1)-C1)*erf((x-q(2))/2/sqrt(10^q(3)))+C1;
            case 'B'
                C = (q(1)-C1)*erf((-x+q(2))/2/sqrt(10^q(3)))+C1;
            case 'C'
                C = (q(1)-C1)*erfc((x-q(2))/2/sqrt(10^q(3)))+C1;
            case 'D'
                C = (q(1)-C1)*erfc((-x+q(2))/2/sqrt(10^q(3)))+C1;
            case 'E'
                C = (q(1)-C1)*(erf((q(4)+x-q(2))/2/sqrt(10^q(3)))+erf((q(4)-x+q(2))/2/sqrt(10^q(3)))-1)+C1;
            case 'F'
                C = (q(1)-C1)*(erfc((q(4)+x-q(2))/2/sqrt(10^q(3)))+erfc((q(4)-x+q(2))/2/sqrt(10^q(3))))+C1;
            case 'G'
                C = (q(1)-C1)/2*(1+erf((x-q(2))/2/sqrt(10^q(3))))+C1;
            case 'H'
                C = (q(1)-C1)/2*erfc((x-q(2))/2/sqrt(10^q(3)))+C1;
        end
    end
% 3. only know C2 ==============================================================
% C1=q(1); x0=q(2); log10(Dt)=q(3); h=q(4);
if isempty(C1) && ~isempty(C2)
    fun = @(q) sqrt(weights).*(f_convolution3(x,model_x,props_matrix,pmsum,profiletype,q,C2)-C);
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'E','F'})
            q0 = [min(C),x0,guess+g,h0];
            lb=[0, min(x), -50, 0];
            ub=[max(C), max(x), 20, inf];
        else
            q0 = [min(C),x0,guess+g];
            lb=[0, min(x), -50];
            ub=[max(C), max(x), 20];
        end
        try
            [q,resnorm,R,~,~,~,J] = lsqnonlin(fun,q0,lb,ub,opts);
            ci = nlparci(q,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=[q(1),C2,q(2:end)];
                ci_min=[ci(1,:);C2,C2;ci(2:end,:)];
                C_pred_min = f_profile3(x,q,profiletype,C2);
                [~,model_y(:,1),model_y(:,2)] = f_convolution3(x,model_x,props_matrix,pmsum,profiletype,q,C2);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution3(x,model_x,props_matrix,pmsum,profiletype,p,C2)
        original=f_profile3(model_x,p,profiletype,C2);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile3(x,q,profiletype,C2)
        switch profiletype
            case 'A'
                C = (C2-q(1))*erf((x-q(2))/2/sqrt(10^q(3)))+q(1);
            case 'B'
                C = (C2-q(1))*erf((-x+q(2))/2/sqrt(10^q(3)))+q(1);
            case 'C'
                C = (C2-q(1))*erfc((x-q(2))/2/sqrt(10^q(3)))+q(1);
            case 'D'
                C = (C2-q(1))*erfc((-x+q(2))/2/sqrt(10^q(3)))+q(1);
            case 'E'
                C = (C2-q(1))*(erf((q(4)+x-q(2))/2/sqrt(10^q(3)))+erf((q(4)-x+q(2))/2/sqrt(10^q(3)))-1)+q(1);
            case 'F'
                C = (C2-q(1))*(erfc((q(4)+x-q(2))/2/sqrt(10^q(3)))+erfc((q(4)-x+q(2))/2/sqrt(10^q(3))))+q(1);
            case 'G'
                C = (C2-q(1))/2*(1+erf((x-q(2))/2/sqrt(10^q(3))))+q(1);
            case 'H'
                C = (C2-q(1))/2*erfc((x-q(2))/2/sqrt(10^q(3)))+q(1);
        end
    end
% 4. known C1 and C2 ==============================================================
% x0=q(1); log10(Dt)=q(2); h=q(3);
if ~isempty(C1) && ~isempty(C2)
    fun = @(q) sqrt(weights).*(f_convolution4(x,model_x,props_matrix,pmsum,profiletype,q,C1,C2)-C);
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'E','F'})
            q0 = [x0,guess+g,h0];
            lb=[min(x), -50, 0];
            ub=[max(x), 20, inf];
        else
            q0 = [x0,guess+g];
            lb=[min(x), -50];
            ub=[max(x), 20];
        end
        try
            [q,resnorm,R,~,~,~,J] = lsqnonlin(fun,q0,lb,ub,opts);
            ci = nlparci(q,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=[C1,C2,q];
                ci_min=[C1,C1;C2,C2;ci];
                C_pred_min = f_profile4(x,q,profiletype,C1,C2);
                [~,model_y(:,1),model_y(:,2)] = f_convolution4(x,model_x,props_matrix,pmsum,profiletype,q,C1,C2);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution4(x,model_x,props_matrix,pmsum,profiletype,p,C1,C2)
        original=f_profile4(model_x,p,profiletype,C1,C2);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile4(x,q,profiletype,C1,C2)
        switch profiletype
            case 'A'
                C = (C2-C1)*erf((x-q(1))/2/sqrt(10^q(2)))+C1;
            case 'B'
                C = (C2-C1)*erf((-x+q(1))/2/sqrt(10^q(2)))+C1;
            case 'C'
                C = (C2-C1)*erfc((x-q(1))/2/sqrt(10^q(2)))+C1;
            case 'D'
                C = (C2-C1)*erfc((-x+q(1))/2/sqrt(10^q(2)))+C1;
            case 'E'
                C = (C2-C1)*(erf((q(3)+x-q(1))/2/sqrt(10^q(2)))+erf((q(3)-x+q(1))/2/sqrt(10^q(2)))-1)+C1;
            case 'F'
                C = (C2-C1)*(erfc((q(3)+x-q(1))/2/sqrt(10^q(2)))+erfc((q(3)-x+q(1))/2/sqrt(10^q(2))))+C1;
            case 'G'
                C = (C2-C1)/2*(1+erf((x-q(1))/2/sqrt(10^q(2))))+C1;
            case 'H'
                C = (C2-C1)/2*erfc((x-q(1))/2/sqrt(10^q(2)))+C1;
        end
    end
end