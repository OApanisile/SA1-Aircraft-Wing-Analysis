function [infa, infb] = panelinf(xa, ya, xb, yb, x, y)
   
    % Vector from A to B
    dx = xb - xa;
    dy = yb - ya;

    L = sqrt(dx^2 + dy^2);  % Panel length
    
    % Unit tangent and normal vectors
    t = [dx, dy] / L;
    n = [-dy, dx] / L;  % Rotate tangent vector 90 degrees CCW
    
    % Vector from A to field point
    rx = x - xa;
    ry = y - ya;

    
    % Transform field point to panel coordinates
    X = rx * t(1) + ry * t(2);
    Yin = rx * n(1) + ry * n(2);

    if abs(Yin) < 1e-5
        Y = 1e-5;
    else
        Y = Yin;
    end

    % Calculating I0
    I0 = -(1 / (4 * pi)) * ( ...
        X*log(X^2 + Y^2) - ...
        (X-L)*log((X-L)^2 + Y^2) - ...
        2*L + ...
        2*Y*(atan(X/Y) - atan((X - L) / Y)) );

    % Calculating I1
    I1 = 1/(8 * pi) * ((X^2 + Y^2) * log(X^2 + Y^2) ...
    - ((X - L)^2 + Y^2) * log((X - L)^2 + Y^2) ...
    - 2 * X * L + L^2);

    % Influence coefficient fa
    infa = (I0*(1-(X/L)) - (I1/L));

    % Influence coefficient fb
    infb = (I0*(X/L) + (I1/L));
  
end