function [int, ils, itr, its, delstar, theta] = bl_solv(x,cp)

    global Re ue0 duedx

    n = length(x); % defines number of panels
    
    ue = sqrt(1-cp);
    duedxov = zeros(size(x));
    duedxov(1) = ue(1)/x(1);
    for i = 2:n
        duedxov(i) = (ue(i) - ue(i-1))/(x(i) - x(i-1));
    end

    f = ueintbit(0, 0, x(1), ue(1));

    m = zeros(size(x));
    H = zeros(size(x));
    He = zeros(size(x));

    theta = zeros(size(x));
    theta(1) = sqrt(f * (0.45/Re) / (ue(1)^6));
    Re_theta = zeros(size(x));
    deltae = zeros(size(x));
    delstar = zeros(size(x));

    int = 0; % Natural transition location
    ils = 0; % Laminar separation location
    itr = 0; % Turbulent reattachment location
    its = 0; % Turbulent separation location
    
    laminar = true; % initializes boundary layer state flag
    i = 1;
    He(1) = 1.57258;
    while laminar && i < n
        i = i + 1;
        
        f = f + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
        
        theta(i) = sqrt(f * (0.45/Re) / (ue(i)^6));
        Re_theta(i) = Re*ue(i)*theta(i);
    
        m(i) = -Re*((theta(i))^2)*(duedxov(i));
        H(i) = thwaites_lookup(m(i));
        He(i) = laminar_He(H(i));

        delstar(i) = H(i) * theta(i);
        
        % Transition test
        if log(Re_theta(i)) >= 18.4*He(i) - 21.74
            laminar = false;
            int = i;
        elseif m(i) >= 0.09
            laminar = false;
            ils = i;
            He(i) = 1.51509;
        end
    end
    
    deltae = He .* theta;
    
    while its == 0 && i < n
        thick0(1) = theta(i);
        thick0(2) = deltae(i);
        ue0 = ue(i);

        i = i + 1;
        duedx = duedxov(i);
    
        [~, thickhist] = ode45(@thickdash,[0, x(i) - x(i-1)],thick0);
    
        theta(i) = thickhist(end,1);
        deltae(i) = thickhist(end,2);
        
        He(i) = deltae(i)/theta(i);
        H(i) = (11*He(i) +15)/(48*He(i) - 59);
        
        
        Re_theta(i) = Re*ue(i)*theta(i);
    
        if ils > 0 && itr == 0 && He(i) > 1.58
    
            itr = i;
    
        elseif He(i) < 1.46
    
            its = i;
            H(i) = 2.803;
    
        end

        delstar(i) = H(i) * theta(i);
    end
    
    He(i:n) = He(i);
    H(i:n) = H(i);
    
    while i < n
        i = i + 1;
        theta(i) = theta(i - 1)*(ue(i - 1)/ue(i)).^(H(i) + 2);
        delstar(i) = H(i) * theta(i);
    end
end

% There is an error in theta such that camber solution is wrong but
% thickness is correct. Error exists from first panel onwards (bl analyser
% works - something deeper at play.