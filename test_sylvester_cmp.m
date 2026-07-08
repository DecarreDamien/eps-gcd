format short 
N = 100; 
d_max = 60;

results_S = NaN(1, N);
time_S = zeros(1, N);
results_G = zeros(1, N);
time_G = zeros(1, N);
degrees = zeros(1, N);

fail_count = 0;

for t = 1:N
    m = 2 + floor(rand(1)*(d_max-1));
    n = 2 + floor(rand(1)*(m-1));
    k = 1 + floor(rand(1)*(n-1));

    deg = m + n;
    degrees(t) = deg;

    p = [1 randn(m,1)'];
    q = [1 randn(n,1)'];

    tic;
    [P, Q, r, f] = epsgcdS(p,q,k);
    time_S(t) = toc;

    if r
        results_S(t) = sup(norm([P-p, Q-q], 2));
    else
        fail_count = fail_count + 1;
    end

    tic;
    [~,~,~,res,~] = uvGCDfixedDegree(p,q,k);
    time_G(t) = toc;
    results_G(t) = res;
end


unique_deg = unique(degrees);

mean_time_S     = zeros(size(unique_deg));
mean_time_G     = zeros(size(unique_deg));
mean_diff       = zeros(size(unique_deg));
std_mean_diff   = zeros(size(unique_deg));

for i = 1:length(unique_deg)
    d = unique_deg(i);
    idx = (degrees == d);

    mean_time_S(i) = mean(time_S(idx), 'omitnan');
    mean_time_G(i) = mean(time_G(idx));

    valid = idx & ~isnan(results_S);
    mean_diff(i) = mean(results_S(valid) - results_G(valid), 'omitnan');
    std_mean_diff(i) = std(results_S(valid) - results_G(valid), 'omitnan');
end

figure;

% Mean time
subplot(2,1,1);
plot(unique_deg, mean_time_S, '-o', 'LineWidth', 1.5); hold on;
plot(unique_deg, mean_time_G, '-s', 'LineWidth', 1.5);
set(gca, 'YScale', 'log');
grid on;
xlabel('Degree (m + n)');
ylabel('Mean time (s)');
title('Mean time vs degré');
legend('epsgcdS', 'uvGCDfixedDegree', 'Location', 'northwest');

% Diff
subplot(2,1,2);
errorbar(unique_deg, mean_diff, std_mean_diff, '-d', 'LineWidth', 1.5);
grid on;
xlabel('Degree (m + n)');
ylabel('mean diff');


fprintf('\nNb failures : %d / %d\n', fail_count, N);