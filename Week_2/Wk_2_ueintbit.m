% Define the ueintbit function
function f = ueintbit(xa, ua, xb, ub)
    ubar = (ua + ub)/2;
    delu = ub - ua;
    delx = xb - xa;

    f = delx*(ubar^5 + (5*(ubar^3)*(delu^2))/6 + (ubar*(delu^4))/16);
end