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