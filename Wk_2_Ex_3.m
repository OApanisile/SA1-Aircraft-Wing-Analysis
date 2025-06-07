xa = 1;
ua = 5;
xb = 3;
ub = 10;
Re_l = 10^5;
dUe_dx = -0.5;
f1 = 0;
f2 = 0;
int = 0; % Natural transition location
ils = 0; % Laminar separation location

x = linspace(0, 1, 101);
ue = 1 * ones(size(x));  % Assuming constant ue = 1
theta = zeros(size(x));
theta_b = zeros(size(x));

% Compute theta using integration
%for i = 2:length(x)
%    
%    f1 = f1 + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
%    theta(i) = sqrt(f1 * (0.45/Re_l) / (ue(i)^6));  % Basic trapezoidal approximation
%
%end

for i = 2:length(x)
    theta_b(i) = (0.664/Re_l^(1/2))*(x(i)^(1/2));
end

ue = uecalc(dUe_dx,ue,x);

n = 101; % defines number of panels
laminar = true; % initializes boundary layer state flag
i = 1;

Re_theta = zeros(1,n);

while laminar && i < n
    i = i + 1;
    f2 = f2 + ueintbit(x(i-1), ue(i-1), x(i), ue(i));
    theta(i) = sqrt(f2 * (0.45/Re_l) / (ue(i)^6));
    Re_theta(i) = Re_l*ue(i)*theta(i);
    m = -Re_l*((theta(i))^2)*(dUe_dx);
    H = thwaites_lookup(m);
    He = laminar_He(H);
   
    if log(Re_theta(i)) >= 18.4*He - 21.74
        laminar = false;
        int = x(i);
        Re_thetant = Re_theta(i);
    
    elseif m >= 0.09
        laminar = false;
        ils = x(i);
        Re_thetals = Re_theta(i);
    end
end

if int ~= 0
    disp(['Natural transition at ' num2str(int) ...
    ' with Re_theta ' num2str(Re_thetant)])
elseif ils ~= 0
    disp(['Laminar separation at ' num2str(ils) ...
    ' with Re_theta ' num2str(Re_thetals)])
end

figure
plot(x, theta, 'b-', x, theta_b, 'r-')
xlabel('x')
ylabel('\theta')
legend('Thwaites', 'Blasius')
title('Momentum thickness \theta vs x')
grid on

% Needs checked - don't like the fact that laminar separation allegedly
% occurs at same x regardless of Re.