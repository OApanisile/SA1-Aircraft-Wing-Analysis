nv = 100;
nx = 51;
ny = 41;
c1 = -1.75:0.25:1.75;
xa = 4.1; 
xb = 2.2; 
ya = 1.3;
yb = 2.9;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
alpha = 0;

% Step (i): Define the cylinder panels
np = 100;
theta = (0:np) * 2 * pi / np;
xs = cos(theta);  % x-coordinates of panel endpoints
ys = sin(theta);  % y-coordinates of panel endpoints

% Step (ii): Define vortex sheet strength at panel edges
g = -2 * sin(theta);  % vortex sheet strength

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

function lhsmat = build_lhs(xs, ys)
    np = length(xs) - 1;             % Number of panels
    psip = zeros(np, np+1);          % Intermediate matrix

    % Loop over control points
    for i = 1:np
        xi = xs(i);                  
        yi = ys(i);                  % Control point at panel edge

        % Loop over all panels
        for j = 1:np
            xj1 = xs(j);
            yj1 = ys(j);
            xj2 = xs(j+1);
            yj2 = ys(j+1);

            % Get influence coefficients
            [infa, infb] = panelinf(xi, yi, xj1, yj1, xj2, yj2);

            % Fill appropriate columns (note: influences from gamma_j and gamma_{j+1})
            psip(i,j)   = psip(i,j)   + infa;  % influence from gamma_j
            psip(i,j+1) = psip(i,j+1) + infb;  % influence from gamma_{j+1}
        end
    end

    % Now add two additional rows for auxiliary conditions
    lhsmat = zeros(np+1, np+1);

    % Copy calculated part (np equations)
    lhsmat(1:np,:) = psip;

    % Final two rows: closure conditions
    % Eq (7): gamma(np+1) = gamma(1)
    lhsmat(np+1,1) = 1;
    lhsmat(np+1,end) = -1;

    % Eq (8): sum of all gamma = 0 (for circulation-free flow)
    % Add if required, or impose different Kutta condition for airfoils
end

lhsmat = build_lhs(xs, ys);



function rhsvec = build_rhs(xs,ys,alpha)
    np = length(xs) - 1;             % Number of panels
    psifs = zeros(np+1, 1);          % Intermediate vector

    % Loop over control points
    for i = 1:np

        % Loop over all panels
        for j = 1:np
            xj1 = xs(j);
            yj1 = ys(j);
            xj2 = xs(j+1);
            yj2 = ys(j+1);

            % Get freestream psi
            psifs1 = yj1*cos(alpha) - xj1*sin(alpha);
            psifs2 = yj2*cos(alpha) - xj2*sin(alpha);

            % Fill vector
            psifs(i, 1) = psifs1 - psifs2;
        end
    end

    % Initialising final vector
    rhsvec = zeros(101, 1);

    % Copy calculated part (np equations)
    rhsvec(1:np+1,1) = psifs;

    % Final two rows: closure conditions
    rhsvec(1,1) = 0;
    rhsvec(end,1) = 0;
end

rhsvec = build_rhs(xs,ys,alpha);

A = build_lhs(xs,ys);
b = build_rhs(xs,ys,alpha);
gam = A\b;