function y = non_causal_filter_v3(b, a, x, alpha)
    poles = roots(a);
    a_stable = poly(poles(abs(poles)<alpha))';
    a_unstable = a(1)*poly(poles(abs(poles)>=alpha))';

    b_unstable = zeros(length(a_unstable),1); b_unstable(end) = 1;
    x1 = flip(filter(b_unstable, flip(a_unstable), flip(x)));
    y = filter(b, a_stable, x1);
end