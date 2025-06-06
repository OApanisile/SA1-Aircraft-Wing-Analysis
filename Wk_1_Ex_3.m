nv = 100;
nx = 51;
ny = 41;
c1 = -0.15:0.05:0.15;
xa = 4.1; 
xb = 2.2; 
ya = 1.3;
yb = 2.9;
xmin = 0;
xmax = 5;
ymin = 0;
ymax = 4;

% Vector from A to B
dx = xb - xa;
dy = yb - ya;

L = sqrt(dx^2 + dy^2);  % Panel length

 
% Generate grid and compute influence cofficients
for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        [infa(i,j), infb(i,j)] = panelinf(xa, ya, xb, yb, xm(i,j), ym(i,j));
    end
end

% Plot fa
figure(1)
contour(xm', ym', infa', c1)
xlabel('x'), ylabel('y')
title('Influence Coefficient f_a')
axis equal

% Plot fb in a separate figure
figure(2)
contour(xm', ym', infb', c1)
xlabel('x'), ylabel('y')
title('Influence Coefficient f_b')
axis equal

yc = 0;
Gamma = 3.0;
c2 = -0.4:0.2:1.2;

function psixy = psipv(xc, yc, Gamma, x, y)
    r_2 = (x - xc)^2 + (y - yc)^2;   % Distance from vortex center
    if r_2 == 0
        error('Field point (x, y) cannot be exactly at the vortex location (xc, yc).');
    end
    psixy = -Gamma / (4 * pi) * log(r_2);
end 

for k = 1:nv
    s = (k - 0.5) / nv;  % Midpoints of segments
    xc(k) = xa + s * (xb - xa);
    yc(k) = ya + s * (yb - ya);
end

% Generate grid and compute streamfunction
for k = 1:nv
    g_a = 1;
    g_b = 0;
    Gamma(k) = (k*(g_a - g_b)/nv)*(L/nv);

    for i = 1:nx
        for j = 1:ny
            psi1(i,j,k) = psipv(xc(nv+1-k), yc(nv+1-k), Gamma(k), xm(i,j), ym(i,j));
        end
    end
end

% Generate grid and compute streamfunction
for k = 1:nv
    g_a = 0;
    g_b = 1;
    Gamma(k) = (k*(g_b - g_a)/nv)*(L/nv);

    for i = 1:nx
        for j = 1:ny
            psi2(i,j,k) = psipv(xc(k), yc(k), Gamma(k), xm(i,j), ym(i,j));
        end
    end
end

infa_approx = sum(psi1, 3);
infb_approx = sum(psi2, 3);

% Plot streamfunction contours
figure(3)
contour(xm', ym', infa_approx', c1)  % Use transpose psi' for correct orientation
xlabel('x'), ylabel('y')
title('Approximate Influence Coefficient f_a')
axis equal

% Plot streamfunction contours
figure(4)
contour(xm', ym', infb_approx', c1)  % Use transpose psi' for correct orientation
xlabel('x'), ylabel('y')
title('Approximate Influence Coefficient f_b')
axis equal