function y=smooth(x, N)
    M = size(x,2);
    y = nan*zeros(size(x));
    
    for i = 1:M
        N2 = 1+round((N-1)/2);
        y(N2:end-N+N2,i) = conv(x(:,i),ones(N,1),'valid')/N;
    end
end