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

% Precompute vortex positions and strengths for fa and fb approximations
g_a = 1;
g_b = 0;
for k = 1:nv
    xc_fa(k) = (del/nv)/2 + (k - 1)*(del/nv);
    Gamma_fa(k) = (k*(g_a - g_b)/nv)*(del/nv);
end

% Compute fa approximation
for k = 1:nv
    for i = 1:nx
        for j = 1:ny
            psi1(i,j,k) = psipv(xc_fa(nv+1-k), yc, Gamma_fa(k), xm(i,j), ym(i,j));
        end
    end
end

g_a = 0;
g_b = 1;
for k = 1:nv
    xc_fb(k) = (del/nv)/2 + (k - 1)*(del/nv);
    Gamma_fb(k) = (k*(g_b - g_a)/nv)*(del/nv);
end

% Compute fb approximation
for k = 1:nv
    for i = 1:nx
        for j = 1:ny
            psi2(i,j,k) = psipv(xc_fb(k), yc, Gamma_fb(k), xm(i,j), ym(i,j));
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