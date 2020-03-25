%% Start
clear;


% Use the modal estimate of beta
beta = [-0.000946480452001148, 0.00295552326330776, 1.30896864696286];


% Reality checks with the data
data = readtable("./data/intermediate/convolution-data.csv");
rows = strcmp("session success rate", data.Metric);
data = data(rows, :);

sampling_variance = mean(data.StdErrorDelta .^ 2);
total_variance = var(data.Delta);
s = beta(2);


% Normal parameters fit with moments
m_moments = mean(data.Delta);
s_moments = sqrt(total_variance - sampling_variance);


%% Marginal distributions of signals
n_points = 20;
delta_grid = linspace(-0.3, 0.3, n_points);

% t distribution
g = @t_distribution;
marginal = @(x) mean(Twee.L(x, data.StdErrorDelta, beta, g));

m = zeros(1, n_points);
for ii = 1:n_points
    m(ii) = marginal(delta_grid(ii));
end

f_t = g(delta_grid, beta);
m_t = m;

close all;
histogram(data.Delta, 'Normalization', 'pdf');
hold on;
plot(delta_grid, f_t);
plot(delta_grid, m);


% Normal with t parameters
g = @(x, beta) normpdf(x, beta(1), beta(2));
marginal = @(x) mean(Twee.L(x, data.StdErrorDelta, beta, g));

m = zeros(1, n_points);
for ii = 1:n_points
    m(ii) = marginal(delta_grid(ii));
end

f_normal_t_param = g(delta_grid, beta);
m_normal_t_param = m;

close all;
histogram(data.Delta, 'Normalization', 'pdf');
hold on;
plot(delta_grid, f_t);
plot(delta_grid, m);


% Normal with moment parameters
g = @(x, beta) normpdf(x, beta(1), beta(2));
marginal = @(x) mean(Twee.L(x, data.StdErrorDelta, [m_moments, s_moments], g));

m = zeros(1, n_points);
for ii = 1:n_points
    m(ii) = marginal(delta_grid(ii));
end

f_normal_t_param = g(delta_grid, [m_moments, s_moments]);
m_normal_t_param = m;

close all;
histogram(data.Delta, 'Normalization', 'pdf');
hold on;
plot(delta_grid, f_t);
plot(delta_grid, m);


% CONCLUSION:
% The normal with t parameters has a nice visual fit. But it fits badly in
% the tails, and for that reason the total variance is less than the
% variance we need to match moments.
% The normal that fits moments has a worse visual fit. The reason is that
% it ends up not having enough mass close to zero. And it still fits the
% data badly, although it evens out these two problems so that it has the
% right amount of variance.
% So it seems appropriate to look at the counterfactual of same parameters
% as the t for the normal distribuiton. Someone making their estimates
% based on robust MLE type things would find similar parameters.


%% Comparing posterior means
n_points = 200;
delta_grid = linspace(-0.15, 0.15, n_points);

mean_sigma = 100;
n = 20e6;
std_error = mean_sigma / sqrt(n);

% Normal
P_normal = @(x) ...
    ((beta(2)^2) .* x + (std_error^2) .* beta(1)) ./ ...
    (beta(2)^2 + std_error^2);

P_normal_moments = @(x) ...
    ((s_moments^2) .* x + (std_error^2) .* m_moments) ./ ...
    (s_moments^2 + std_error^2);

% t
g = @t_distribution;
P_t = @(x) ...
    Twee.mean_posterior(x, std_error, beta, g);

y_normal = zeros(n_points, 1);
y_normal_moments = zeros(n_points, 1);
y_t = zeros(n_points, 1);
for ii = 1:n_points
    delta = delta_grid(ii);
    y_normal(ii) = P_normal(delta);
    y_normal_moments(ii) = P_normal_moments(delta);
    y_t(ii) = P_t(delta);
end

% graph
font_size = 16;
close all;
figure();
hold on;
plot(delta_grid, y_t, 'k', 'LineWidth', 2);
plot(delta_grid, y_normal, 'k--', 'LineWidth', 2);
plot(delta_grid, y_normal_moments, 'k.', 'LineWidth', 2);
legend({...
    'Student t', ...
    'Normal', ...
    'Normal (matched moments)'}, ...
    'Location', 'northwest', 'FontSize', font_size);
legend('boxoff');
xlabel('Estimated quality $\hat \delta _i$ in the experiment', 'Interpreter' , 'Latex');
ylabel('Posterior mean quality $P(\hat \delta _i)$', 'Interpreter' , 'Latex');
ax = gca;
ax.FontSize = font_size;


%% Thresholds
% Zero
t_star_zero = fzero(P_t, 0) / std_error
p_value_zero = 1 - normcdf(t_star_zero)

% Positive implementation cost
c = 0.05;
f = @(x) P_t(x) - c;
t_star_c = fzero(f, 0) / std_error
p_value_c = 1 - normcdf(t_star_c)
