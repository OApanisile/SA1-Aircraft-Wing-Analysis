xc = 0.75;
yc = 0.5;
Gamma = 3.0;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
nx = 51;
ny = 41;
c = -0.4:0.2:1.2;

function psixy = psipv(xc, yc, Gamma, x, y)
    r_2 = (x - xc)^2 + (y - yc)^2;   % Distance from vortex center
    if r_2 == 0
        error('Field point (x, y) cannot be exactly at the vortex location (xc, yc).');
    end
    psixy = -Gamma / (4 * pi) * log(r_2);
end

% Generate grid and compute streamfunction
for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        psi(i,j) = psipv(xc, yc, Gamma, xm(i,j), ym(i,j));
    end
end

% Plot streamfunction contours
contour(xm, ym, psi, c)  % Use transpose psi' for correct orientation
xlabel('x'), ylabel('y')
title('Streamfunction of a Point Vortex')
axis equal