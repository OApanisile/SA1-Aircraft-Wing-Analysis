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
 
% Generate grid and compute influence cofficients
for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        psi(i,j) = ym(i,j);
        for p = 1:np
            [infa(i,j), infb(i,j)] = panelinf(xs(p), ys(p), xs(p+1), ys(p+1), xm(i,j), ym(i,j));
            psi(i,j) = psi(i,j) + (g(p)*infa(i,j)) + (g(p+1)*infb(i,j));
        end
    end
end

% Plot psi
figure(1)
contour(xm', ym', psi', c1)
xlabel('x'), ylabel('y')
title('Influence Coefficient Psi')
axis equal
