function y = mean1(x, Mx)
    N = size(x,1);
    M = size(x,2);

    xv = zeros(Mx,1);
    y = zeros(N,M);
    for m=1:M
        for n=1:N
            xv = [x(n,m); xv(1:end-1)];
            y(n,m) = mean(xv,1);
        end
    end
end