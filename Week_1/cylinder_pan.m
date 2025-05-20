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

% Step (i): Define the cylinder panels
np = 100;
theta = (0:np) * 2 * pi / np;
xs = cos(theta);  % x-coordinates of panel endpoints
ys = sin(theta);  % y-coordinates of panel endpoints

% Step (ii): Define vortex sheet strength at panel edges
g = -2 * sin(theta);  % vortex sheet strength

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

    [infa,infb] = refpaninf(L,X,Y);
  
end

function lhsmat = build_lhs(xs,ys)
    %build_lhs Obtains the lhs of the streamfunction and gamma equation
    %   
    
    np = length(xs) - 1;
    psip = zeros(np,np+1); %Initialise matrix to reduce computation times
    
    for i=1:np
        for j=1:np+1    
                
                if j == 1 
                    [infa,infb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
                    psip(i,j) = infa;
                elseif j == np+1 %Accounts for index out of bounds
                    [infa1,infb1] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
                    psip(i,j) = infb1;
                else
                    [infa,infb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
                    [infa1,infb1] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
                    psip(i,j) = infa + infb1;
                end
        end
    end %Looping over all the points around the cylinder measuring the contribution from each`point vortex onto one single point of the cylinder
    
    %Kutta Condition
    
    lhsmat = zeros(np+1,np+1);
    lhsmat(1, [1:3,np-1:np]) = [2, -2, 1, -1, 2];
    lhsmat(np+1,[2:3,np-1:np+1]) = [-2, 1, -1, 2, -2];
    
    
    for i=2:np
        for j=1:np+1
            lhsmat(i,j) = psip(i,j) - psip(i-1,j);
        end
    end %Set our lhs matrix to the difference of the streamfunctions contribution due to the pannels of two consecutive points
end

function rhsvec = build_rhs(xs,ys,alpha) 
    %build_rhs Obtains the rhs the Streamfunction and gamma equation
    %   Detailed explanation goes here
    
    np = length(xs) - 1; %Number of points we are dividing the cylinder into
    rhsvec = zeros(np+1,1); %Initialise the vector to reduce computation times
    psi = zeros(np+1,1);
    
    for i = 1:np+1
        psi(i) = ys(i)*cos(alpha) - xs(i)*sin(alpha); %Free stream contribution to the streamfunction
    end
    
    for i = 2:np
        rhsvec(i) = psi(i-1) - psi(i); %The rhs of our equation will be the difference between the free stream contributions between two consecutive points
    end

end

lhsmat = build_lhs(xs,ys);
rhsvec = build_rhs(xs,ys,alpha); %Unsure 0s in right place

A = build_lhs(xs,ys);
b_0 = build_rhs(xs,ys,0);
b_a = build_rhs(xs,ys,alpha);

gam_0 = A\b_0;
gam_a = A\b_a;

circulation_0 = 0;
circulation_a = 0;

L1 = 2*pi/(np);
x = theta(2);
L = sqrt(2*(1-cos(x)));
for i = 1:np+1
    circulation_0 = circulation_0 + gam_0(i)*L;
    circulation_a = circulation_a + gam_a(i)*L;
end

% Plot psi
figure(1)
plot(theta/pi,gam_0)
xlabel('x'), ylabel('y')
title('Influence Coefficient Psi')
axis equal

% Plot psi
figure(2)
plot(theta/pi,gam_a)
xlabel('x'), ylabel('y')
title('Influence Coefficient Psi')
axis equal

for i = 1:np
    for j = 1:np
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(np-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(np-1);
    end
end

% Plot psi
figure(3)
contour(xm', ym', A', c1)
xlabel('x'), ylabel('y')
title('Influence Coefficient Psi')
axis equal