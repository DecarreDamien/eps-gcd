format short 
N = 500;
d_max = 50;
depth_max = 8;

succes_rate = 0;
for i=1:N
    m = 2 + floor(rand(1)*(d_max-1));
    n = 2 + floor(rand(1)*(m-1));
    k = 1 + floor(rand(1)*(n-1));
    g = [1 randn(k,1)'];
    p = conv(g, [1 randn(m-k+1,1)']);
    q = conv(g, [1 randn(n-k+1,1)']);
    e_max = 8;
    e_min = -12;
    e = uncertainty(e_min, e_max);
    %e = midrad(0, 10^(-2));
    [r, E] = epsgcdB(p+e, q+e, k, 10);
    for j=1:depth_max
        if r 
            break 
        else 
            if E.code == 1
                e_max = (e_max + e_min) / 2;
            elseif E.code == 2
                e_min = (e_max + e_min) / 2;
            else 
                break 
            end
            e = uncertainty(e_min, e_max);
            %e = e/10;
            [r, E] = epsgcdB(p+e, q+e, k, 10);
        end
    end
    if r 
        succes_rate = succes_rate + 1;
    end
end
succes_rate = succes_rate / N

function x = uncertainty(k1, k2)
    x = midrad(0, 10^((k1 + k2)/2));
end