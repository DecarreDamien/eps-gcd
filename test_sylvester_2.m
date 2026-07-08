function [norm_dis, degrees, ex, flag, failure] ...
    = test_sylvester_error(d_max, N, delta)

    failure = 0;
    norm_dis = zeros(1, N);
    degrees = zeros(1, N);
    ex = zeros(1, N);
    flag = zeros(1, N);
    for i=1:N
        m = 2 + floor(rand(1)*(d_max-1));
        n = 2 + floor(rand(1)*(m-1));
        k = 1 + floor(rand(1)*(n-1));
        degrees(i) = m + n;
        g = [1 randn(k,1)'];
        p = conv(g, [1 randn(m-k,1)']) + delta;
        q = conv(g, [1 randn(n-k,1)']) - delta;
        [P, Q, r, f] = epsgcdS(p,q,k);
        flag(i) = f;
        if r 
            ex(i) = 1;
        end
        if ~isnan(P(1))
            norm_valp = sup(norm(P-p, inf));
            norm_valq = sup(norm(Q-q, inf));
            norm_dis(i) = maxi(norm_valp, norm_valq);
        else 
            failure = failure + 1;
        end
    end 
end

function [m] = maxi(a,b)
    if a > b 
        m = a;
    else 
        m = b;
    end
end


N = 1000; d_max = 50;
[dis10, deg10, ex10, flag10, fr10] = test_sylvester_error(d_max, N, 1e-10);
fr10 = fr10 / N
[dis6, deg6, ex6, flag6, fr6] = test_sylvester_error(d_max, N, 1e-6);
fr6 = fr6 / N
[dis2, deg2, ex2, flag2, fr2] = test_sylvester_error(d_max, N, 1e-2);
fr2 = fr2 / N


%% Radius Norm distribution when k is the deg of the gcd of f and g
%subplot(1,3,1);
%bins = 10.^(floor(log10(min(dis10(ex10 == 1), [], "all"))): ...
%    ceil(log10(max(dis10(ex10 == 1), [], "all"))));
%h1 = histogram(dis10(ex10 == 1), bins);
%set(gca, "XScale", "log");
%ylabel('Occurrence');
%xlabel('Output radius with an initial perturbation of 10^{-10}');
%subplot(1,3,2);
%bins = 10.^(floor(log10(min(dis6(ex6 == 1), [], "all"))): ...
%    ceil(log10(max(dis6(ex6 == 1), [], "all"))));
%h2 = histogram(dis6(ex6 == 1), bins);
%set(gca, "XScale", "log");
%xlabel('Output radius with an initial perturbation of 10^{-6}');
%subplot(1,3,3);
%bins = 10.^(floor(log10(min(dis2(ex2 == 1), [], "all"))): ...
%    ceil(log10(max(dis2(ex2 == 1), [], "all"))));
%h3 = histogram(dis2(ex2 == 1), bins);
%set(gca, "XScale", "log");
%xlabel('Output radius with an initial perturbation of 10^{-2}');
