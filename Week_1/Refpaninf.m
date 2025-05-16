del = 1.5;
X = 2;
Yin = 1.5;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
nx = 51;
ny = 41;
c = -0.15:0.05:0.15;
function [infa, infb] = refpaninf(del, X, Yin)
    
    if abs(Yin) < 1e-5
        Y = 1e-5;
    else
        Y = Yin;
    end

    % Precompute reused terms
    X1 = X;
    X2 = X - del;
    Y2 = Y^2;

    % Influence coefficient fa
    infa = -(1 / (4 * pi)) * ( ...
        X*log(X^2 + Y^2) - ...
        (X-del)*log((X-del)^2 + Y^2) - ...
        2*del + ...
        2*Y*(atan(X/Y) - atan((del - X) / Y)) );

    % Influence coefficient fb
    infb = 1/(8 * pi) * ((X^2 + Y^2) * log(X^2 + Y^2) ...
    - ((X - del)^2 + Y^2) * log((X - del)^2 + Y^2) ...
    - 2 * X * del + del^2);
  
end
 
% Generate grid and compute influence cofficients
for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        [infa(i,j), infb(i,j)] = refpaninf(del, xm(i,j), ym(i,j));
    end
end

% Plot fa
figure(1)
contour(xm', ym', infa', c)
xlabel('x'), ylabel('y')
title('Influence Coefficient f_a')
axis equal

% Plot fb in a separate figure
figure(2)
contour(xm', ym', infb', c)
xlabel('x'), ylabel('y')
title('Influence Coefficient f_b')
axis equal