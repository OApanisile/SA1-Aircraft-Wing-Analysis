del = 1.5;
X = 2;
Yin = 1.5;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
nv = 100;
nx = 51;
ny = 41;
c1 = -0.15:0.05:0.15;

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
    Y = rx * n(1) + ry * n(2);

    
    if abs(Yin) < 1e-5
        Y = 1e-5;
    else
        Y = Yin;
    end

    % Precompute reused terms
    X1 = X;
    X2 = X - del;
    Y2 = Y^2;

    % Calculating I0
    I0 = -(1 / (4 * pi)) * ( ...
        X*log(X^2 + Y^2) - ...
        (X-del)*log((X-del)^2 + Y^2) - ...
        2*del + ...
        2*Y*(atan(X/Y) - atan((X - del) / Y)) );

    % Calculating I1
    I1 = 1/(8 * pi) * ((X^2 + Y^2) * log(X^2 + Y^2) ...
    - ((X - del)^2 + Y^2) * log((X - del)^2 + Y^2) ...
    - 2 * X * del + del^2);

    % Influence coefficient fa
    infa = (I0*(1-(X/del)) - (I1/del));

    % Influence coefficient fb
    infb = (I0*(X/del) + (I1/del));
  
end
% Panel endpoints
xa = 4.1; ya = 1.3;
xb = 2.2; yb = 2.9;

% Domain grid for field points
[x, y] = meshgrid(linspace(1, 6, 100), linspace(0, 5, 100));
fa = zeros(size(x));
fb = zeros(size(x));

% Evaluate influence at each point
for i = 1:size(x,1)
    for j = 1:size(x,2)
        [fa(i,j), fb(i,j)] = panelinf(xa, ya, xb, yb, x(i,j), y(i,j));
    end
end

% Plot exact contours
figure;
subplot(2,2,1);
contour(x, y, fa, 50);
title('Exact f_a');
axis equal; colorbar;

subplot(2,2,2);
contour(x, y, fb, 50);
title('Exact f_b');
axis equal; colorbar;

% Approximate using multiple point vortices (e.g., 10)
n = 10;
s = linspace(0, 1, n);
xp = xa + (xb - xa) * s;
yp = ya + (yb - ya) * s;

% Assume uniform strength distribution and compute total influence
fap = zeros(size(x));
fbp = zeros(size(x));

for k = 1:n
    rx = x - xp(k);
    ry = y - yp(k);
    r2 = rx.^2 + ry.^2;
    fa_k = -ry ./ (2*pi*r2);
    fb_k = rx ./ (2*pi*r2);
    fap = fap + fa_k;
    fbp = fbp + fb_k;
end

subplot(2,2,3);
contour(x, y, fap, 50);
title('Approx. f_a (point vortices)');
axis equal; colorbar;

subplot(2,2,4);
contour(x, y, fbp, 50);
title('Approx. f_b (point vortices)');
axis equal; colorbar;
