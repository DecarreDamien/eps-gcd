format short
m = 4; n = 3; k = 2;
g = [1 randn(k,1)'];
p = conv(g, [1 randn(m-k,1)'])
q = conv(g, [1 randn(n-k,1)'])

%default version
%p = p + midrad(0, e)
%q = q + midrad(0, e)
e = 1e-8;
%adding small perturbation so that p and q are not in the
%center of the interval polynomials
alpha_p = e * rand(1, m+1);
beta_p  = e * rand(1, m+1);
alpha_q = e * rand(1, n+1);
beta_q  = e * rand(1, n+1);

P = infsup(p - alpha_p, p + beta_p)
Q = infsup(q - alpha_q, q + beta_q)

[~, E] = epsgcdB(P, Q, k, 10)