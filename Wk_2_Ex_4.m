clear

global ue0 duedx Re;

Re = 10^7;
duedx = 0;
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
plot(x, theta7, 'b-', x, theta9, 'r-', x, theta, 'g-')
xlabel('x')
ylabel('\theta')
legend('1/7th Law', '1/9th Law', 'Calculated')
title('Momentum thickness \theta vs x')
grid on