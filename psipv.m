function psixy = psipv(xc, yc, Gamma, x, y)
    r_2 = (x - xc)^2 + (y - yc)^2;   % Distance from vortex center
    if r_2 == 0
        error('Field point (x, y) cannot be exactly at the vortex location (xc, yc).');
    end
    psixy = -Gamma / (4 * pi) * log(r_2);
end
