format short 
e_values = [midrad(0,1e-14), midrad(0,1e-10), midrad(0,1e-6),...
    midrad(0,1e-2)];

results = zeros(5, length(e_values));

for ie = 1:length(e_values)
    e = e_values(ie);
    res = zeros(5,1);

    for k = 2:21
        for m = k+1:k+15
            for n = k+1:k+15
                g = poly(randn(k,1));
                p = conv(g, poly(randn(m-k+1,1)));
                q = conv(g, poly(randn(n-k+1,1)));

                [~, E] = epsgcdB(p+e, q+e, k, 10);
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

labels = {'Succès', 'B_{m-k} not full rank', 'M_i not full rank',...
    'no suitable X', 'B_{m-k+1} full rank'};

xticks(1:length(labels));
xticklabels(labels);

legend({'e=1e-14','e=1e-10','e=1e-6','e=1e-2'});

