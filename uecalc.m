function ue = uecalc(dUe_dx,ue,x)
    ue(1) = 1;
    for n = 2:length(x)
        ue(n) = ue(n-1) + (dUe_dx*0.01);
    end
end 
