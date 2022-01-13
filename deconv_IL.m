function [p_min, ci_min, C_pred_min, model_y] = deconv_IL(x, C, C1, C2, C3, profiletype, weights, factor, model_x, props_matrix, pmsum)
% Set up fittype and options.
% C1=p(1); C2=p(2); x0=p(3); log10(Dt)=p(4); h=p(5); C3=p(6);
p_min=[]; ci_min=[]; C_pred_min=[]; model_y=[]; rsn_min=1e10;
%initial guess
h0=(max(x)-min(x))/2;
x0=(max(x)+min(x))/2;
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
% 1. unknown C1, C2 and C3==============================================================
if isempty(C1) && isempty(C2) && isempty(C3)
    fun = @(p) sqrt(weights).*(f_convolution1(x,model_x,props_matrix,pmsum,profiletype,p)-C);
    if strcmp(profiletype,'K')
        C3=max(C(1),C(end)); %initial C3 guess
    elseif strcmp(profiletype,'L')
        C3=min(C(1),C(end)); %initial C3 guess
    end
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'I','J'})
            p0 = [min(C),max(C),x0,guess+g,h0];
            lb=[0, min(C), min(x), -50, 0];
            ub=[max(C), inf, max(x), 20, max(x)-min(x)];
        else
            p0 = [min(C),max(C),x0,guess+g,h0,C3];
            lb=[0, min(C), min(x), -50, 0, 0];
            ub=[max(C), inf, max(x), 20, max(x)-min(x), max(C)];
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
            case 'I'
                C = (p(2)-p(1))/2*(erf((p(5)-x+p(3))/2/sqrt(10^p(4)))+erf((p(5)+x-p(3))/2/sqrt(10^p(4))))+p(1);
            case 'J'
                C = (p(2)-p(1))/2*(erfc((p(5)-x+p(3))/2/sqrt(10^p(4)))+erfc((p(5)+x-p(3))/2/sqrt(10^p(4))))+p(1);
            case 'K'
                C = (p(2)-p(6))/2*erf((p(5)-x+p(3))/2/sqrt(10^p(4)))+(p(2)-p(1))/2*erf((p(5)+x-p(3))/2/sqrt(10^p(4)))+(p(1)+p(6))/2;
            case 'L'
                C = (p(6)-p(1))/2*erfc((p(5)-x+p(3))/2/sqrt(10^p(4)))+(p(2)-p(1))/2*erfc((p(5)+x-p(3))/2/sqrt(10^p(4)))+p(1);
        end
    end
% 2. only know C1 ==============================================================
% C2=q(1); x0=q(2); log10(Dt)=q(3); h=q(4); C3=q(5);
if ~isempty(C1) && isempty(C2) && isempty(C3)
    fun = @(q) sqrt(weights).*(f_convolution2(x,model_x,props_matrix,pmsum,profiletype,q,C1)-C);
    if strcmp(profiletype,'K')
        C3=max(C(1),C(end)); %initial C3 guess
    elseif strcmp(profiletype,'L')
        C3=min(C(1),C(end)); %initial C3 guess
    end
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'I','J'})
            q0 = [max(C),x0,guess+g,h0];
            lb=[min(C), min(x), -50, 0];
            ub=[inf, max(x), 20, max(x)-min(x)];
        else
            q0 = [max(C),x0,guess+g,h0,C3];
            lb=[min(C), min(x), -50, 0, 0];
            ub=[inf, max(x), 20, max(x)-min(x), max(C)];
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
            case 'I'
                C = (q(1)-C1)/2*(erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+erf((q(4)+x-q(2))/2/sqrt(10^q(3))))+C1;
            case 'J'
                C = (q(1)-C1)/2*(erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+erfc((q(4)+x-q(2))/2/sqrt(10^q(3))))+C1;
            case 'K'
                C = (q(1)-q(5))/2*erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+(q(1)-C1)/2*erf((q(4)+x-q(2))/2/sqrt(10^q(3)))+(C1+q(5))/2;
            case 'L'
                C = (q(5)-C1)/2*erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+(q(1)-C1)/2*erfc((q(4)+x-q(2))/2/sqrt(10^q(3)))+C1;
        end
    end
% 3. only know C2 ==============================================================
% C1=q(1); x0=q(2); log10(Dt)=q(3); h=q(4); C3=q(5);
if isempty(C1) && ~isempty(C2) && isempty(C3)
    fun = @(q) sqrt(weights).*(f_convolution3(x,model_x,props_matrix,pmsum,profiletype,q,C2)-C);
    if strcmp(profiletype,'K')
        C3=max(C(1),C(end)); %initial C3 guess
    elseif strcmp(profiletype,'L')
        C3=min(C(1),C(end)); %initial C3 guess
    end
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'I','J'})
            q0 = [min(C),x0,guess+g,h0];
            lb=[0, min(x), -50, 0];
            ub=[max(C), max(x), 20, max(x)-min(x)];
        else
            q0 = [min(C),x0,guess+g,h0,C3];
            lb=[0, min(x), -50, 0, 0];
            ub=[max(C), max(x), 20, max(x)-min(x), max(C)];
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
            case 'I'
                C = (C2-q(1))/2*(erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+erf((q(4)+x-q(2))/2/sqrt(10^q(3))))+q(1);
            case 'J'
                C = (C2-q(1))/2*(erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+erfc((q(4)+x-q(2))/2/sqrt(10^q(3))))+q(1);
            case 'K'
                C = (C2-q(5))/2*erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+(C2-q(1))/2*erf((q(4)+x-q(2))/2/sqrt(10^q(3)))+(q(1)+q(5))/2;
            case 'L'
                C = (q(5)-q(1))/2*erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+(C2-q(1))/2*erfc((q(4)+x-q(2))/2/sqrt(10^q(3)))+q(1);
        end
    end
% 4. only know C1 and C2 ==============================================================
% x0=q(1); log10(Dt)=q(2); h=q(3); C3=q(4);
if ~isempty(C1) && ~isempty(C2) && isempty(C3)
    fun = @(q) sqrt(weights).*(f_convolution4(x,model_x,props_matrix,pmsum,profiletype,q,C1,C2)-C);
    if strcmp(profiletype,'K')
        C3=max(C(1),C(end)); %initial C3 guess
    elseif strcmp(profiletype,'L')
        C3=min(C(1),C(end)); %initial C3 guess
    end
    % Fit model to data.
    for g=gg
        if contains(profiletype, {'I','J'})
            q0 = [x0,guess+g,h0];
            lb=[min(x), -50, 0];
            ub=[max(x), 20, max(x)-min(x)];
        else
            q0 = [x0,guess+g,h0,C3];
            lb=[min(x), -50, 0, 0];
            ub=[max(x), 20, max(x)-min(x), max(C)];
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
            case 'I'
                C = (C2-C1)/2*(erf((q(3)-x+q(1))/2/sqrt(10^q(2)))+erf((q(3)+x-q(1))/2/sqrt(10^q(2))))+C1;
            case 'J'
                C = (C2-C1)/2*(erfc((q(3)-x+q(1))/2/sqrt(10^q(2)))+erfc((q(3)+x-q(1))/2/sqrt(10^q(2))))+C1;
            case 'K'
                C = (C2-q(4))/2*erf((q(3)-x+q(1))/2/sqrt(10^q(2)))+(C2-C1)/2*erf((q(3)+x-q(1))/2/sqrt(10^q(2)))+(C1+q(4))/2;
            case 'L'
                C = (q(4)-C1)/2*erfc((q(3)-x+q(1))/2/sqrt(10^q(2)))+(C2-C1)/2*erfc((q(3)+x-q(1))/2/sqrt(10^q(2)))+C1;
        end
    end
% 5. only know C3 ==============================================================
% C1=q(1); C2=q(2); x0=q(3); log10(Dt)=q(4); h=q(5);
if isempty(C1) && isempty(C2) && ~isempty(C3)
    fun = @(q) sqrt(weights).*(f_convolution5(x,model_x,props_matrix,pmsum,profiletype,q,C3)-C);
    % Fit model to data.
    for g=gg
        q0 = [min(C),max(C),x0,guess+g,h0];
        lb=[0, min(C), min(x), -50, 0];
        ub=[max(C), inf, max(x), 20, max(x)-min(x)];
        try
            [q,resnorm,R,~,~,~,J] = lsqnonlin(fun,q0,lb,ub,opts);
            ci = nlparci(q,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=[q,C3];
                ci_min=[ci;C3,C3];
                C_pred_min = f_profile5(x,q,profiletype,C3);
                [~,model_y(:,1),model_y(:,2)] = f_convolution5(x,model_x,props_matrix,pmsum,profiletype,q,C3);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution5(x,model_x,props_matrix,pmsum,profiletype,p,C3)
        original=f_profile5(model_x,p,profiletype,C3);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile5(x,q,profiletype,C3)
        switch profiletype
            case 'I'
                C = (q(2)-q(1))/2*(erf((q(5)-x+q(3))/2/sqrt(10^q(4)))+erf((q(5)+x-q(3))/2/sqrt(10^q(4))))+q(1);
            case 'J'
                C = (q(2)-q(1))/2*(erfc((q(5)-x+q(3))/2/sqrt(10^q(4)))+erfc((q(5)+x-q(3))/2/sqrt(10^q(4))))+q(1);
            case 'K'
                C = (q(2)-C3)/2*erf((q(5)-x+q(3))/2/sqrt(10^q(4)))+(q(2)-q(1))/2*erf((q(5)+x-q(3))/2/sqrt(10^q(4)))+(q(1)+C3)/2;
            case 'L'
                C = (C3-q(1))/2*erfc((q(5)-x+q(3))/2/sqrt(10^q(4)))+(q(2)-q(1))/2*erfc((q(5)+x-q(3))/2/sqrt(10^q(4)))+q(1);
        end
    end
% 6. only know C1 and C3 ==============================================================
% C2=q(1); x0=q(2); log10(Dt)=q(3); h=q(4);
if ~isempty(C1) && isempty(C2) && ~isempty(C3)
    fun = @(q) sqrt(weights).*(f_convolution6(x,model_x,props_matrix,pmsum,profiletype,q,C1,C3)-C);
    % Fit model to data.
    for g=gg
            q0 = [max(C),x0,guess+g,h0];
            lb=[min(C), min(x), -50, 0];
            ub=[inf, max(x), 20, max(x)-min(x)];
        try
            [q,resnorm,R,~,~,~,J] = lsqnonlin(fun,q0,lb,ub,opts);
            ci = nlparci(q,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=[C1,q,C3];
                ci_min=[C1,C1;ci;C3,C3];
                C_pred_min = f_profile6(x,q,profiletype,C1,C3);
                [~,model_y(:,1),model_y(:,2)] = f_convolution6(x,model_x,props_matrix,pmsum,profiletype,q,C1,C3);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution6(x,model_x,props_matrix,pmsum,profiletype,p,C1,C3)
        original=f_profile6(model_x,p,profiletype,C1,C3);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile6(x,q,profiletype,C1,C3)
        switch profiletype
            case 'I'
                C = (q(1)-C1)/2*(erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+erf((q(4)+x-q(2))/2/sqrt(10^q(3))))+C1;
            case 'J'
                C = (q(1)-C1)/2*(erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+erfc((q(4)+x-q(2))/2/sqrt(10^q(3))))+C1;
            case 'K'
                C = (q(1)-C3)/2*erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+(q(1)-C1)/2*erf((q(4)+x-q(2))/2/sqrt(10^q(3)))+(C1+C3)/2;
            case 'L'
                C = (C3-C1)/2*erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+(q(1)-C1)/2*erfc((q(4)+x-q(2))/2/sqrt(10^q(3)))+C1;
        end
    end
% 7. only know C2 and C3 ==============================================================
% C1=q(1); x0=q(2); log10(Dt)=q(3); h=q(4);
if isempty(C1) && ~isempty(C2) && ~isempty(C3)
    fun = @(q) sqrt(weights).*(f_convolution7(x,model_x,props_matrix,pmsum,profiletype,q,C2,C3)-C);
    % Fit model to data.
    for g=gg
            q0 = [min(C),x0,guess+g,h0];
            lb=[0, min(x), -50, 0];
            ub=[max(C), max(x), 20, max(x)-min(x)];
        try
            [q,resnorm,R,~,~,~,J] = lsqnonlin(fun,q0,lb,ub,opts);
            ci = nlparci(q,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=[q(1),C2,q(2:end),C3];
                ci_min=[ci(1,:);C2,C2;ci(2:end,:);C3,C3];
                C_pred_min = f_profile7(x,q,profiletype,C2,C3);
                [~,model_y(:,1),model_y(:,2)] = f_convolution7(x,model_x,props_matrix,pmsum,profiletype,q,C2,C3);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution7(x,model_x,props_matrix,pmsum,profiletype,p,C2,C3)
        original=f_profile7(model_x,p,profiletype,C2,C3);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile7(x,q,profiletype,C2,C3)
        switch profiletype
            case 'I'
                C = (C2-q(1))/2*(erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+erf((q(4)+x-q(2))/2/sqrt(10^q(3))))+q(1);
            case 'J'
                C = (C2-q(1))/2*(erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+erfc((q(4)+x-q(2))/2/sqrt(10^q(3))))+q(1);
            case 'K'
                C = (C2-C3)/2*erf((q(4)-x+q(2))/2/sqrt(10^q(3)))+(C2-q(1))/2*erf((q(4)+x-q(2))/2/sqrt(10^q(3)))+(q(1)+C3)/2;
            case 'L'
                C = (C3-q(1))/2*erfc((q(4)-x+q(2))/2/sqrt(10^q(3)))+(C2-q(1))/2*erfc((q(4)+x-q(2))/2/sqrt(10^q(3)))+q(1);
        end
    end
% 8. known C1, C2, and C3 ==============================================================
% x0=q(1); log10(Dt)=q(2); h=q(3);
if ~isempty(C1) && ~isempty(C2) && ~isempty(C3)
    fun = @(q) sqrt(weights).*(f_convolution8(x,model_x,props_matrix,pmsum,profiletype,q,C1,C2,C3)-C);
    % Fit model to data.
    for g=gg
            q0 = [x0,guess+g,h0];
            lb=[min(x), -50, 0];
            ub=[max(x), 20, max(x)-min(x)];
        try
            [q,resnorm,R,~,~,~,J] = lsqnonlin(fun,q0,lb,ub,opts);
            ci = nlparci(q,R,'jacobian',J,'alpha',0.05); %parameter bounds at 95c.l.
            if resnorm<rsn_min
                rsn_min=resnorm;
                p_min=[C1,C2,q,C3];
                ci_min=[C1,C1;C2,C2;ci;C3,C3];
                C_pred_min = f_profile8(x,q,profiletype,C1,C2,C3);
                [~,model_y(:,1),model_y(:,2)] = f_convolution8(x,model_x,props_matrix,pmsum,profiletype,q,C1,C2,C3);
            end
        catch
        end
    end
end
    function [convoluted_downsampled,convoluted,original]=f_convolution8(x,model_x,props_matrix,pmsum,profiletype,p,C1,C2,C3)
        original=f_profile8(model_x,p,profiletype,C1,C2,C3);
        convoluted=original.*props_matrix';
        convoluted=sum(convoluted');
        convoluted=convoluted.*pmsum;
        convoluted=convoluted(((size(convoluted,2)-size(model_x,2))/2)+1:size(convoluted,2)-((size(convoluted,2)-size(model_x,2))/2));
        convoluted_downsampled=interp1(model_x,convoluted,x,'nearest');
    end
    function C =f_profile8(x,q,profiletype,C1,C2,C3)
        switch profiletype
            case 'I'
                C = (C2-C1)/2*(erf((q(3)-x+q(1))/2/sqrt(10^q(2)))+erf((q(3)+x-q(1))/2/sqrt(10^q(2))))+C1;
            case 'J'
                C = (C2-C1)/2*(erfc((q(3)-x+q(1))/2/sqrt(10^q(2)))+erfc((q(3)+x-q(1))/2/sqrt(10^q(2))))+C1;
            case 'K'
                C = (C2-C3)/2*erf((q(3)-x+q(1))/2/sqrt(10^q(2)))+(C2-C1)/2*erf((q(3)+x-q(1))/2/sqrt(10^q(2)))+(C1+C3)/2;
            case 'L'
                C = (C3-C1)/2*erfc((q(3)-x+q(1))/2/sqrt(10^q(2)))+(C2-C1)/2*erfc((q(3)+x-q(1))/2/sqrt(10^q(2)))+C1;
        end
    end

end