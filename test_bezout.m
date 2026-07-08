format short 
e_values = [1e-14, 1e-10, 1e-6, 1e-2];

results = zeros(5, length(e_values));

for ie = 1:length(e_values)
    e = e_values(ie);
    res = zeros(5,1);

    for k = 2:21
        for m = k+1:k+15
            for n = k+1:k+15
                g = [1 randn(k,1)'];
                p = conv(g, [1 randn(m-k,1)']);
                q = conv(g, [1 randn(n-k,1)']);
                
                %default version
                %P = p + midrad(0, e)
                %Q = q + midrad(0, e)
                
                %adding small perturbation so that p and q are not in the
                %center of the interval polynomials
                alpha_p = e * rand(1, m+1);
                beta_p  = e * rand(1, m+1);
                alpha_q = e * rand(1, n+1);
                beta_q  = e * rand(1, n+1);
                
                P = infsup(p - alpha_p, p + beta_p);
                Q = infsup(q - alpha_q, q + beta_q);

                [~, E] = epsgcdB(P, Q, k, 10);
                i = E.code + 1;
                res(i) = res(i) + 1;
            end
        end
    end

    results(:, ie) = res;
end


figure;
bar(results, 'grouped');
colormap(parula(length(e_values)));

labels = {'Success', 'B_{m-k} not full rank', 'M_i not full rank',...
    'no suitable X', 'B_{m-k+1} full rank'};

xticks(1:length(labels));
xticklabels(labels);

legend({'e=1e-14','e=1e-10','e=1e-6','e=1e-2'});

