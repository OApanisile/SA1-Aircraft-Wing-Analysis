xa = 1;
ua = 5;
xb = 3;
ub = 10;
Re_l = 2500;
dUe_dx = 0.2;

x = linspace(0, 1, 101);
ue = 1 * ones(size(x));  % Assuming constant ue = 1
theta = zeros(size(x));
theta_b = zeros(size(x));

% Define the ueintbit function
function f = ueintbit(xa, ua, xb, ub)
    ubar = (ua + ub)/2;
    delu = ub - ua;
    delx = xb - xa;

    f = delx*(ubar^5 + (5*(ubar^3)*(delu^2))/6 + (ubar*(delu^4))/16);
end

f = 0;

% Compute theta using integration
for i = 2:length(x)
    xa = x(i-1);
    xb = x(i);
    ua = ue(i-1);
    ub = ue(i);
    
    f = f + ueintbit(xa, ua, xb, ub);
    theta(i) = sqrt(f * (0.45/Re_l) / (ue(i)^6));  % Basic trapezoidal approximation
end

for i = 2:length(x)
    theta_b(i) = (0.664/Re_l^(1/2))*(x(i)^(1/2));
end

function ue = uecalc(dUe_dx,ue,x)
    ue(1) = 1;
    for n = 2:length(x)
        ue(n) = ue(n-1) + (dUe_dx*0.01); %X spacing.
    end
end 

ue = uecalc(dUe_dx,ue,x);

f = 0;

% Compute theta using integration
for i = 2:length(x)
    xa = x(i-1);
    xb = x(i);
    ua = ue(i-1);
    ub = ue(i);
    
    f = f + ueintbit(xa, ua, xb, ub);
    theta_2(i) = sqrt(f * (0.45/Re_l) / (ue(i)^6));  % Basic trapezoidal approximation
end

for i = 2:length(x)
    theta_b_2(i) = (0.664/Re_l^(1/2))*(x(i)^(1/2));
end

Re_theta = zeros(size(x));
m = zeros(size(x));

for i = 1:length(x)
    Re_theta(i) = Re_l*ue(i)*theta(i);
    m(i) = -1*Re_l*(theta(i)^2)*(dUe_dx);
end 

function H = thwaites_lookup(m)
%
%  function H = thwaites_lookup(m)
%
%  Returns interpolated value of shape factor 
%  for the specified value of Thwaites' 'm' parameter.  Based on tabulated
%  values in Young, p90.
%

if m <= -0.25     % outside range of table, return values at -0.25

%  l = 0.5;
  H = 2;

elseif m <= 0.09    % interpolate from tabulated values

  mtab = [-.25 -.2 -.14 -.12 -.1 -.08 -.064 -.048 -.032 -.016 .0 .016 .032 ...
             .04 .048 .056 .06 .064 .068 .072 .076 .080 .084 .086 .088 .09];

%  ltab = [.5 .463 .404 .382 .359 .333 .313 .291 .268 .244 .22 .195 .168 ... 
%            .153 .138 .122 .113 .104 .095 .085 .072 .056 .038 .027 .015 .0];

  Htab = [2 2.07 2.18 2.23 2.28 2.34 2.39 2.44 2.49 2.55 2.61 2.67 2.75 ...
           2.81 2.87 2.94 2.99 3.04 3.09 3.15 3.22 3.3 3.39 3.44 3.49 3.55];

%  l = interp1(mtab,ltab,m);
  H = interp1(mtab,Htab,m);

else      %  m > 0.09; b-l separated, return m = 0.09 values

%  l = 0;
  H = 3.55;

end
end

function He = laminar_He(H)
%
%  function He = laminar_He(H)
%
%  Provides energy shape factor He as function of shape factor H.
%  Based on Eppler & Somers' expressions for H(He).
%  For H < 1.855:              He = 1.742
%      1.855 < H < 2.591       1.742 > He > 1.573 (analytical inversion)
%      2.591 < H < 4.029       1.573 > He > 1.515 (iterative inversion)
%

if H <= 1.855 

  He = 89.582142/(2*25.715786);

elseif H <= 2.5911

  He = (89.582142 - sqrt(89.582142^2 - 4*25.715786*(79.870845-H))) ...
                                                          /(2*25.715786);

else

  He = 1.545;
  Heold = 0;

  while abs(He-Heold) > 0.0005

    dum = (4.02922 - H)/(583.60182 - 724.55916*He + 227.1822*He^2);
    Heold = He;
    He = 1.51509 + dum^2;

  end

end
end


He = laminar_He(thwaites_lookup(m));

figure(1)
plot(x, theta, 'b-', x, theta_b, 'r-')
xlabel('x')
ylabel('\theta')
legend('Thwaites', 'Blasius')
title('Momentum thickness \theta vs x')
grid on

figure(2)
plot(x, theta_2, 'b-', x, theta_b_2, 'r-')
xlabel('x')
ylabel('\theta')
legend('Thwaites', 'Blasius')
title('Momentum thickness \theta vs x')
grid on