del = 1.5;
xmin = -2.5;
xmax = 2.5;
ymin = -2;
ymax = 2;
nv = 100;
nx = 51;
ny = 41;
c = -0.15:0.05:0.15;
yc = 0;

function [infa, infb] = refpaninf(del, X, Yin)
    
    if abs(Yin) < 1e-5
        Y = 1e-5;
    else
        Y = Yin;
    end

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

function psixy = psipv(xc, yc, Gamma, x, y)
    r_2 = (x - xc)^2 + (y - yc)^2;   % Distance from vortex center
    if r_2 == 0
        error('Field point (x, y) cannot be exactly at the vortex location (xc, yc).');
    end
    psixy = -Gamma / (4 * pi) * log(r_2);
end 

% Generate grid and compute streamfunction
for k = 1:nv
    xc(k) = (del/nv)/2 + (k - 1)*(del/nv);
    g_a = 1;
    g_b = 0;
    Gamma(k) = (k*(g_a - g_b)/nv)*(del/nv);
    for i = 1:nx
        for j = 1:ny
            psi1(i,j,k) = psipv(xc(nv+1-k), yc, Gamma(k), xm(i,j), ym(i,j));
        end
    end
end


for k = 1:nv
    xc(k) = (del/nv)/2 + (k - 1)*(del/nv);
    g_a = 0;
    g_b = 1;
    Gamma(k) = (k*(g_b - g_a)/nv)*(del/nv);
    for i = 1:nx
        for j = 1:ny
            psi2(i,j,k) = psipv(xc(k), yc, Gamma(k), xm(i,j), ym(i,j));
        end
    end
end

infa_approx = sum(psi1, 3);
infb_approx = sum(psi2, 3);

% Plot streamfunction contours
figure(3)
contour(xm', ym', infa_approx', c)  % Use transpose psi' for correct orientation
xlabel('x'), ylabel('y')
title('Approximate Influence Coefficient f_a')
axis equal

% Plot streamfunction contours
figure(4)
contour(xm', ym', infb_approx', c)  % Use transpose psi' for correct orientation
xlabel('x'), ylabel('y')
title('Approximate Influence Coefficient f_b')
axis equal