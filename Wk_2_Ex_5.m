clear

global Re ue0 duedx;

Re = 10^8;
duedx = -0.6;
ue0 = 1;

x0 = 0.01;
thick0(1) = 0.037*x0*(Re*x0)^(-1/5);
thick0(2) = 1.80*thick0(1);

[delx, thickhist] = ode45(@thickdash,[0 0.99],thick0);

theta = thickhist(:,1);
deltae = thickhist(:,2);

x = x0 + delx;

theta7 = 0.037*x.*(Re*x).^(-1/5);
theta9 = 0.023*x.*(Re*x).^(-1/6);

figure(1)
plot(x, theta, 'b-', x, deltae, 'r-')
xlabel('x')
ylabel('\theta')
legend('Momentum thickness', 'Energy thickness')
title('Momentum thickness \theta vs x') % Not sure what to title plot
grid on

He = deltae./theta;
He_sep = 1.46*ones(1,length(x));

figure (2)
plot(x, He, 'b-', x, He_sep, 'r-')
xlabel('x')
ylabel('He')
legend('Energy shape factor He', 'Turbulent boundary layer separation threshold')
title('Energy shape factor He vs x')
grid on