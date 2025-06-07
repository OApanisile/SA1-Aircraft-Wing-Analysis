xa = 1;
ua = 5;
xb = 3;
ub = 10;
Re_l = 100*10^6;
dUe_dx = 0.2;
f1 = 0;
f2 = 0;

x = linspace(0, 1, 101);
ue = 1 * ones(size(x));  % Assuming constant ue = 1
theta = zeros(size(x));
theta_b = zeros(size(x));
theta_2 = zeros(size(x));

% Compute theta using integration
for i = 2:length(x)
    f1 = f1 + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
    theta(i) = sqrt(f1 * (0.45/Re_l) / (ue(i)^6));
end

for i = 2:length(x)
    theta_b(i) = (0.664/Re_l^(1/2))*(x(i)^(1/2));
end

ue = uecalc(dUe_dx,ue,x);

n = 101; % defines number of panels
laminar = true; % initializes boundary layer state flag
i = 1;
while laminar && i < n
    i = i + 1;
    
    f2 = f2 + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
    theta_2(i) = sqrt(f2 * (0.45/Re_l) / (ue(i)^6));
    Re_theta(i) = Re_l*ue(i)*theta(i);

    m = -Re_l*((theta(i))^2)*(dUe_dx);
    H = thwaites_lookup(m);
    He = laminar_He(H);
    
    if log(Re_theta(i)) >= 18.4*He - 21.74
        laminar = false;
        %xt = x(i);
        %Re_thetat = Re_theta(i)/1000;
        fprintf('For i = %.4f, x(i) = %.4f and Re_theta/1000 = %.4f\n', i, x(i), Re_theta(i)/1000)
        disp([x(i) Re_theta(i)/1000])
    end
end

figure(1)
plot(x, theta, 'b-', x, theta_b, 'r-')
xlabel('x')
ylabel('\theta')
legend('Thwaites', 'Blasius')
title('Momentum thickness \theta vs x')
grid on

figure(2)
plot(x, theta_2, 'b-', x, theta_b, 'r-')
xlabel('x')
ylabel('\theta')
legend('Thwaites', 'Blasius')
title('Momentum thickness \theta vs x')
grid on