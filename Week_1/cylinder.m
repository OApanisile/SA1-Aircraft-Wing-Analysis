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
 
% Generate grid and compute influence cofficients
for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        psi(i,j) = ym(i,j);
        for p = 1:np
            [infa(i,j), infb(i,j)] = panelinf(xs(p), ys(p), xs(p-1), ys(p-1), xm(i,j), ym(i,j));
            psi(i,j) = psi(i,j) + (g(p)*infa(i,j)) + (g(p-1)*infb(i,j));
        end
    end
end

% Plot psi
figure(1)
contour(xm', ym', psi', c1)
xlabel('x'), ylabel('y')
title('Influence Coefficient Psi')
axis equal
