function p = calc_random_path(N)
global min_path_gain;
path_ok = false;
while ~path_ok
    d1 = rand*N/4; d2 = rand*N/4; d3 = rand*N/4;
    p  = (2*rand-1)*sinc((0:N-1)'-d1)+(2*rand-1)*sinc((0:N-1)'-d2)+(2*rand-1)*sinc((0:N-1)'-d3);
    pf = freqz(p, 1, 1024);
    if min(abs(pf)) > min_path_gain
        path_ok = true;
    end
end
