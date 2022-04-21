function [a,b]=adjust_plant(a,b,stability_margin)

        roots_a = roots(a);
        ind = abs(roots_a)>1-stability_margin;
        roots_a(ind) = roots_a(ind)./abs(roots_a(ind))*(1-stability_margin);
        
        roots_b = roots(b);
        ind = abs(roots_b)>1-stability_margin & abs(roots_b)<1;
        roots_b(ind) = roots_b(ind)./abs(roots_b(ind))*(1-stability_margin);
        ind = abs(roots_b)<1+stability_margin & abs(roots_b)>=1;
        roots_b(ind) = roots_b(ind)./abs(roots_b(ind))*(1+stability_margin);
        
        a = poly(roots_a)';
        b = b(1)*poly(roots_b)';
        
        b = b*sqrt(sum(a.^2)/sum(b.^2));

end