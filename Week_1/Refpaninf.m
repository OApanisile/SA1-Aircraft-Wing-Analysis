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

% Generate grid and compute streamfunction
for k = 1:nv
    xc(k) = (del/nv)/2 + (k - 1)*(del/nv);
    for i = 1:nx
        for j = 1:ny
            psi(i,j,k) = psipv(xc(k), yc, Gamma, xm(i,j), ym(i,j));
        end
    end
end

psitot = sum(psi, 3);

psi(:,:,1)

% Plot streamfunction contours
figure(3)
contour(xm', ym', psitot', c2)  % Use transpose psi' for correct orientation
xlabel('x'), ylabel('y')
title('Streamfunction of a Point Vortex')
axis equal