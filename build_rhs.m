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