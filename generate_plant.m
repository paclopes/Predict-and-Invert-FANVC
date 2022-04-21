function [a,b] = generate_plant(Nx, stability_margin)
    a_ok = false;
    while ~a_ok
        a = 1/Nx*randn(Nx+1,1); a(1) = 1;
        if max(abs(roots(a)))<1-stability_margin
            a_ok = true;
        end
    end
    b_ok = false;
    while ~b_ok
        b = randn(Nx+1,1);
        if min(abs(abs(roots(b))-1))>stability_margin
            b_ok = true;
        end
    end
    
    b = b*sqrt(sum(a.^2)/sum(b.^2));
end