xa = 1;
ua = 5;
xb = 3;
ub = 10;
Re_l = 2500;

x = linspace(0, 1, 101);
ue = 1 * ones(size(x));  % Assuming constant ue = 1
theta = zeros(size(x));
theta_b = zeros(size(x));

% Define the ueintbit function
function f = ueintbit(xa, ua, xb, ub)
    ubar = (ua + ub)/2;
    delu = ub - ua;
    delx = xb - xa;

    f = delx*(ubar^5 + (5*(ubar^3)*(delu^2))/6 + (ubar*(delu^4))/16);
end

f = 0;

% Compute theta using integration
for i = 2:length(x)
    xa = x(i-1);
    xb = x(i);
    ua = ue(i-1);
    ub = ue(i);
    
    f = f + ueintbit(xa, ua, xb, ub);
    theta(i) = sqrt(f * (0.45/Re_l) / (ue(i)^6));  % Basic trapezoidal approximation
end

for i = 2:length(x)
    theta_b(i) = (0.664/Re_l^(1/2))*(x(i)^(1/2));
end

figure(1)
plot(x, theta, 'b-', x, theta_b, 'r-')
xlabel('x')
ylabel('\theta')
legend('Thwaites', 'Blasius')
title('Momentum thickness \theta vs x')
grid on