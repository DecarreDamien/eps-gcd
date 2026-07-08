format short 

%m = 32; n = 21; k = 4;

%g = [1 randn(k,1)'];
%p = conv(g, [1 randn(m-k,1)'])+1e-6;
%q = conv(g, [1 randn(n-k,1)'])-1e-6;

%tic
%[u,v,w,res,cond] = uvGCD(p,q,1e-5);
%toc
%disp(res)
%disp(length(u)-1)

%d1 = norm(conv(u,v) - p, inf)
%d2 = norm(conv(u,w) - q, inf)
%P = intval(p) + midrad(0, res);
%all(in(conv(u,v)', P))
%Q = intval(q) + midrad(0, res);
%all(in(conv(u,w)', Q))
%[~, E] = epsgcdB(P,Q,length(u)-1, 10)

N = 100; d_max = 60;

results = zeros(5,1);
degrees = zeros(1, N);
for t=1:N
    m = 2 + floor(rand(1)*(d_max-1));
    n = 2 + floor(rand(1)*(m-1));
    k = 1 + floor(rand(1)*(n-1));
    degrees(t) = m + n;
    g = [1 randn(k,1)'];
    p = conv(g, [1 randn(m-k,1)'])+1e-6;
    q = conv(g, [1 randn(n-k,1)'])-1e-6;

    [u,v,w,res,cond] = uvGCD(p,q,1e-4,0);

    P = intval(p) + midrad(0, res);
    Q = intval(q) + midrad(0, res);
    [~, E] = epsgcdB(P,Q,length(u), 10);
    i = E.code + 1;
    results(i) = results(i) + 1;
end

figure;
bar(results);

labels = {'Success', 'B_{m-k} not full rank', 'M_i not full rank',...
    'no suitable X', 'B_{m-k+1} full rank'};

xticks(1:length(labels));
xticklabels(labels);