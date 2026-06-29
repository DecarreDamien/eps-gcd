function [norm_dis, norm_dis_flag, norm_dis_c, degrees, ex, failure] ...
    = test_sylvester_error(d_max, N)

    failure = 0;
    norm_dis = zeros(1, N);
    norm_dis_flag = zeros(1, N);
    norm_dis_c = zeros(1, N);
    degrees = zeros(1, N);
    ex = zeros(1, N);
    for i=1:N
        m = 2 + floor(rand(1)*(d_max-1));
        n = 2 + floor(rand(1)*(m-1));
        k = 1 + floor(rand(1)*(n-1));
        degrees(i) = m + n;
        g = poly(randn(k,1));
        p = conv(g,poly(randn(m-k+1,1)));
        q = conv(g,poly(randn(n-k+1,1)));
        [P, Q, r, f] = epsgcdS(p,q,k);
        if r 
            ex(i) = 1;
        end
        if ~isnan(P(1))
            norm_valp = sup(norm(P-p, inf));
            norm_valq = sup(norm(Q-q, inf));
            if f >= 0
                norm_dis_flag(i) = maxi(norm_valp, norm_valq);

                [P, Q, ~, ~] = epsgcdS(p,q,f);
                norm_valp = sup(norm(P-p, inf));
                norm_valq = sup(norm(Q-q, inf));
                norm_dis_c(i) = maxi(norm_valp, norm_valq);
            else
                norm_dis(i) = maxi(norm_valp, norm_valq);
            end
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


N = 10000; d_max = 50;
[dis, dis_f, dis_c, deg, ex, fr] = test_sylvester_error(d_max, N);
%fr = fr / N

%% Radius Norm distribution when k is ok
%subplot(1,3,1);
%bins = 10.^(floor(log10(min(dis(dis > 0), [], "all"))): ...
%    ceil(log10(max(dis, [], "all"))));
%h1 = histogram(dis, bins);
%set(gca, "XScale", "log");

%% Radius Norm distribution when k is not ok
%subplot(1,3,2);
%bins_f = 10.^(floor(log10(min(dis_f(dis_f > 0), [], "all"))): ...
%    ceil(log10(max(dis_f, [], "all"))));
%h2 = histogram(dis_f, bins_f);
%set(gca, "XScale", "log");

%% Corrected Radius Norm distribution when k is not ok
%subplot(1,3,3);
%bins_c = 10.^(floor(log10(min(dis_c(dis_c > 0), [], "all"))): ...
%    ceil(log10(max(dis_c, [], "all"))));
%h3 = histogram(dis_c, bins_c);
%set(gca, "XScale", "log");

%% Number of bad k vs degree
%T = table(deg(:), dis(:), 'VariableNames', {'Categorie', 'Valeurs'});
%bilan = groupsummary(T, 'Categorie', @(x) mean(x == 0), 'Valeurs');
%figure;
%bilan.Properties.VariableNames = {'Categorie', 'NombreElements', 'ProportionZeros'};
%bar(bilan.Categorie, bilan.ProportionZeros);

%% Log of Geometric Mean radius vs degree (deg p + deg q)
%T = table(deg(:), dis(:), 'VariableNames', {'Categorie', 'Valeurs'});
%bilan = groupsummary(T, 'Categorie', @(x) mean(log10(x(x > 0))), 'Valeurs');
%figure;
%bilan.Properties.VariableNames = {'Categorie', 'NombreElements', 'Moyenne'};
%bar(bilan.Categorie, bilan.Moyenne);

%% Number of (f,g) s.t. deg gcd(f,g) = k
%T = table(deg(:), ex(:), 'VariableNames', {'Categorie', 'Valeurs'});
%bilan = groupsummary(T, 'Categorie', @(x) mean(x), 'Valeurs');
%figure;
%bilan.Properties.VariableNames = {'Categorie', 'NombreElements', 'ProportionExact'};
%bar(bilan.Categorie, bilan.ProportionExact);