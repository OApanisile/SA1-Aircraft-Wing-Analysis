function lhsmat = build_lhs(xs,ys)
    np = length(xs) - 1;
    psip = zeros(np,np+1);
    
    for i=1:np
        for j=1:np+1    
                if j == 1 
                    [infa,infb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
                    psip(i,j) = infa;
                elseif j == np+1
                    [infa1,infb1] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
                    psip(i,j) = infb1;
                else
                    [infa,infb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
                    [infa1,infb1] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
                    psip(i,j) = infa + infb1;
                end
        end
    end %Looping over all the points around the cylinder measuring the contribution from each point vortex to a single point on the cylinder
    
    %Kutta Condition
    lhsmat = zeros(np+1,np+1);

    % Original Kutta condition:
    lhsmat(1,1) = 1;
    lhsmat(np+1,np+1) = 1;
    
    %lhsmat(1, [1:3,np-1:np]) = [1, -1, 0.5, -0.5, 1];  
    %lhsmat(np+1, [2:3,np-1:np+1]) = [1, -0.5, 0.5, -1, 1];  
    
    
    for i=2:np
        for j=1:np+1
            lhsmat(i,j) = psip(i,j) - psip(i-1,j);
        end
    end %Set our lhs matrix to the difference of the streamfunctions contribution due to the pannels of two consecutive points
end