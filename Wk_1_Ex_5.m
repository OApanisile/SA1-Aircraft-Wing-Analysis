% Parameters
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
alpha = 0.1;

np = 100;
theta = (0:np) * 2 * pi / np;
xs = cos(theta); % x-coordinates of panel endpoints
ys = sin(theta); % y-coordinates of panel endpoints

% --- Build System ---
A = build_lhs(xs, ys);
b_0 = build_rhs(xs, ys, 0);
b_a = build_rhs(xs, ys, alpha);

gam_0 = A \ b_0;
gam_a = A \ b_a;

% --- Compute Circulations ---
circulation_0 = sum(gam_0) * (2*pi/np);
circulation_a = sum(gam_a) * (2*pi/np);

fprintf('circulation_0 = %.4f\n', circulation_0)
fprintf('circulation_a = %.4f\n', circulation_a)


% --- Plot Gamma Distributions ---
figure(1)
plot(theta/pi, gam_0)
xlabel('theta/\pi'), ylabel('\gamma_0')
title('Surface Velocity \gamma vs \theta/\pi for \alpha = 0')
axis([0 2 -2.5 2.5])

figure(2)
plot(theta/pi, gam_a)
xlabel('theta/\pi'), ylabel('\gamma_a')
title('Surface Velocity \gamma vs \theta/\pi for \alpha = 0.1')
axis([0 2 -2.5 2.5])

% --- Compute and Plot Streamfunctions ---
[xg, yg] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny));
psitot_0 = zeros(size(xg));
psitot_a = zeros(size(xg));

% alpha = 0
for i = 1:size(xg,1)
    for j = 1:size(xg,2)
        x = xg(i,j);
        y = yg(i,j);
        if sqrt(x^2 + y^2) < 1
            psitot_0(i,j) = NaN;
            continue
        end
        psi_free = y;
        psi_ind = 0;
        for k = 1:np
            [fa, fb] = panelinf(xs(k), ys(k), xs(k+1), ys(k+1), x, y);
            psi_ind = psi_ind + gam_0(k)*fa + gam_0(k+1)*fb;
        end
        psitot_0(i,j) = psi_free + psi_ind;
    end
end

% alpha = 0.1
for i = 1:size(xg,1)
    for j = 1:size(xg,2)
        x = xg(i,j);
        y = yg(i,j);
        if sqrt(x^2 + y^2) < 1
            psitot_a(i,j) = NaN;
            continue
        end
        psi_free = y*cos(alpha) - x*sin(alpha);
        psi_ind = 0;
        for k = 1:np
            [fa, fb] = panelinf(xs(k), ys(k), xs(k+1), ys(k+1), x, y);
            psi_ind = psi_ind + gam_a(k)*fa + gam_a(k+1)*fb;
        end
        psitot_a(i,j) = psi_free + psi_ind;
    end
end

% Plot alpha = 0
figure(3)
contour(xg, yg, psitot_0, 50)
hold on
plot(cos(theta), sin(theta), 'k')
xlabel('x'), ylabel('y')
title('Total Streamfunction \psi + \psi_p (\alpha = 0)')
axis equal
grid on

% Plot alpha = 0.1
figure(4)
contour(xg, yg, psitot_a, 50)
hold on
plot(cos(theta), sin(theta), 'k')
xlabel('x'), ylabel('y')
title('Total Streamfunction \psi + \psi_p (\alpha = 0.1)')
axis equal
grid on