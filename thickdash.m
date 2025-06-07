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