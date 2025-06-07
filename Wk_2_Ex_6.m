clear

global Re ue0 duedx;
Re = 10^5;
duedx = -0.39; % Critical velocity gradient is between -0.38 and -0.39
ue0 = 1;

f = 0;
int = 0; % Natural transition location
ils = 0; % Laminar separation location
itr = 0; % Turbulent reattachment location
its = 0; % Turbulent separation location

x = linspace(0, 1, 101);
ue = 1 * ones(size(x));  % Assuming constant ue = 1
theta = zeros(size(x));
theta_b = zeros(size(x));

for i = 2:length(x)
    theta_b(i) = (0.664/Re^(1/2))*(x(i)^(1/2));
end

function ue = uecalc(dUe_dx,ue,x)
    ue(1) = 1;
    for n = 2:length(x)
        ue(n) = ue(n-1) + (dUe_dx*0.01); %X spacing.
    end
end 

ue = uecalc(duedx,ue,x);

n = 101; % defines number of panels
laminar = true; % initializes boundary layer state flag
i = 1;
He(1) = 1.57258;
while laminar && i < n
    i = i + 1;
    
    f = f + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
    theta(i) = sqrt(f * (0.45/Re) / (ue(i)^6));
    Re_theta(i) = Re*ue(i)*theta(i);

    m(i) = -Re*((theta(i))^2)*(duedx);
    H(i) = thwaites_lookup(m(i));
    He(i) = laminar_He(H(i));
    
    % Transition test
    if log(Re_theta(i)) >= 18.4*He(i) - 21.74
        laminar = false;

        int = x(i);
        Re_thetant = Re_theta(i);
    elseif m(i) >= 0.09
        laminar = false;

        ils = x(i);
        Re_thetals = Re_theta(i);
        He(i) = 1.51509;
    end
end

deltae(i) = He(i) .* theta(i);

ue0 = ue0 + i/n*duedx;

while its == 0 && i < n
    thick0(1) = theta(i);
    thick0(2) = deltae(i);
    i = i + 1;

    [delx, thickhist] = ode45(@thickdash,[0, x(i) - x(i-1)],thick0);

    ue0 = ue0 +1/n*duedx;
    theta(i) = thickhist(end,1);
    deltae(i) = thickhist(end,2);
    
    He(i) = deltae(i)/theta(i);

    if He(i) >= 1.46
        H(i) = (11*He(i) +15)/(48*He(i) - 59);
    else
        H(i) = 2.803;
    end
    
    Re_theta(i) = Re*ue(i)*theta(i);

    if ils > 0 && itr == 0 && He(i) > 1.58

        itr = x(i);
        Re_thetatr = Re_theta(i);

    elseif He(i) < 1.46

        its = x(i);
        Re_thetats = Re_theta(i);

    end
end

He(i:n) = He(i);
H(i:n) = H(i);

while i < n
    i = i + 1;
    theta(i) = theta(i - 1)*(ue(i - 1)/ue(i)).^(H(i) + 2);
end

if int ~= 0
    disp(['Natural transition at ' num2str(int) ...
    ' with Re_theta ' num2str(Re_thetant)])
end

if ils ~= 0
    disp(['Laminar separation at ' num2str(ils) ...
    ' with Re_theta ' num2str(Re_thetals)])
end

if itr ~= 0
    disp(['Turbulent reattachment at ' num2str(itr) ...
    ' with Re_theta ' num2str(Re_thetatr)])
end

if its ~= 0
    disp(['Turbulent separation at ' num2str(its) ...
    ' with Re_theta ' num2str(Re_thetats)])
end

figure(1)
plot(x, theta, 'b-', x, theta_b, 'r-')
xlabel('x')
ylabel('\theta')
legend('Thwaites', 'Blasius')
title('Momentum thickness \theta vs x')
grid on

figure(2)
plot(x, He, 'r-')
xlabel('x')
ylabel('He')
%legend('Thwaites')
title('Energy shape factor He vs x')
grid on