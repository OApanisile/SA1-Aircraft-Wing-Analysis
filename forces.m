
function [cl, cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)

    % Calculate cl
    cl = -2*circ;

    % Calculate cd
    u_te = sqrt(1-cp(end));
    theta_te = thetal(end) + thetau(end);
    delta_te = delstarl(end) + delstaru(end);
    H_te = delta_te/theta_te;

    theta_inf = theta_te*(u_te)^((H_te+5)/2);
    cd = 2*theta_inf;
end