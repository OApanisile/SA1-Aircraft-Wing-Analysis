clear

global Re ue0 duedx;
Re = 10^8;
duedx = -0.6;
ue0 = 1;

x0 = 0.01;
thick0(1) = 0.037*x0*(Re*x0)^(-1/5);
thick0(2) = 1.80*thick0(1);

function dthickdx = thickdash(xmx0,thick)
    
    % Passing global variables into function
    global Re ue0 duedx

    % Calculate ue
    ue = ue0 + duedx*xmx0;

    % Calculate Re_theta
    Re_theta = Re*ue*thick(1);

    % Calculate He
    He = thick(2)/thick(1);

    % Determining H 
    if He >= 1.46
        H = (11*He +15)/(48*He - 59);
    else
        H = 2.803;
    end

    % Calculate cf
    cf = 0.091448*((H - 1)*Re_theta)^(-0.232) * exp(-1.26*H);

    % Calculate cdiss
    cdiss = 0.010025*((H - 1)*Re_theta)^(-1/6);

    % Set output vector
    dthickdx = [cf/2 - (H + 2)/ue(1) * duedx * thick(1); cdiss - 3/ue *duedx * thick(2)];

end

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