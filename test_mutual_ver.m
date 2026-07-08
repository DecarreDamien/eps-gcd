format short
function [sr] = test_mutual_verification(d_max, N)
    sr = 0;
    for i=1:N
        m = 2 + floor(rand(1)*(d_max-1));
        n = 2 + floor(rand(1)*(m-1));
        k = 1 + floor(rand(1)*(n-1));
        g = [1 randn(k,1)'];
        p = conv(g, [1 randn(m-k+1,1)']);
        q = conv(g, [1 randn(n-k+1,1)']);
        [P, Q, ~, ~] = epsgcdS(p,q,k);
        if epsgcdB(P, Q, k, 100)
            sr = sr + 1;
        end
    end 
end

[succes_rate] = test_mutual_verification(30, 100)
%ouput : 27