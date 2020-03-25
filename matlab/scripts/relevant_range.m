% Additional results requested by JPE reviewers in reject and resubmit
% This file calculates over what ranges the small n and large n
% approximations are good.
%% Start
clear;
close all;
addpath './matlab/classes';
addpath './matlab/functions';


%% Load data
data       = readtable('./data/intermediate/convolution-data.csv');
mle_output = load('./matlab/mat/mle.mat');
mle_output = mle_output.mle_output;
mle_output = mle_output(strcmp(mle_output.algorithm, 'gradient knitro'), :);


%% Create variables
% Prior
g          = @t_distribution;
dg         = @d_t_distribution;

g_normal = @(xx, bb) normpdf(xx, bb(1), bb(2));

% Select data for success rate
data = data(strcmp(data.Metric, 'session success rate'), :);

% Select MLE output with preferred specification of no weights.
row = (mle_output.use_weights == 0);
mle_output = mle_output(row, :);

% Choose beta for success rate with benchmark estimate
row = strcmp(mle_output.metric, 'session success rate');
beta = mle_output(row, 'beta');
beta = table2array(beta);

% Create useful variables
delta_hat = data.Delta;
sigma_n = data.StdErrorDelta;
sigma = data.sigma;

mean_sigma = mean(sigma);

% Sort variables by the signal
[delta_hat, sorted_index] = sort(delta_hat);
sigma = sigma(sorted_index);
sigma_n = sigma_n(sorted_index);

% Set normal distribution betas
beta_normal_shape = beta(1:2);
beta_normal_moments = [0, 0];
beta_normal_moments(1) = mean(delta_hat);
beta_normal_moments(2) = var(delta_hat) - mean(sigma_n .^ 2);
beta_normal_moments(2) = sqrt(beta_normal_moments(2));


%% Plot the production function and the theoretical approximations
% Define variables for the approximations
f_infinity = Twee.f(Inf, mean_sigma, beta, g);
g_zero = g(0, beta);
c_t_infinity = ...
    gamma((beta(3) + 1)/2) / gamma(beta(3)/2) / sqrt(beta(3) * pi) * ...
    (beta(3) * beta(2)^2 )^ ( (beta(3)+1) / 2) / ...
    beta(2); % Power law constant c for the t distribution for large values.
t_one = sqrt(2 * (beta(3) - 1) * log(mean_sigma / sqrt(1))); % For the small n formula, we will assume that the t statistic is constant at a value equal to the approximation at n = 1.

approximation_small_constant = ...
    1 / 2 * beta(3) * c_t_infinity * ...
    (mean_sigma * t_one) ^ (1 - beta(3));

approximation_small = @(n) ...
    2 / (beta(3) - 1) * approximation_small_constant * ...
    n .^ ((beta(3) - 1)/2);

approximation_large = @(n) ...
    f_infinity - ...
    (1/2) * g_zero * mean_sigma ^2 ./ n;


% Define grid for n
% Higher resolution for the empirically relevant values.
% And a logged grid covering a broad swath of values.
n_points_target = 100;
n_grid_small = linspace(1, 1e8, n_points_target);
n_grid_large = linspace(log(1), log(1e11), n_points_target);
n_grid_large = exp(n_grid_large);
n_grid = [n_grid_small, n_grid_large];
n_grid = sort(unique(n_grid));
n_points = length(n_grid);


% Loop for calculations
f_grid = zeros(n_points, 1);
approximation_small_grid = f_grid;
approximation_large_grid = f_grid;
ratio_grid = f_grid; % This is the ratio of f and n^alpha-1/2.
error_small_grid = f_grid;
error_large_grid = f_grid;

for ii = 1:n_points
    n = n_grid(ii);
    f_grid(ii) = Twee.f(n, mean_sigma, beta, g);
    approximation_small_grid(ii) = approximation_small(n);
    approximation_large_grid(ii) = approximation_large(n);
    % approximation_large_grid(ii) = max(approximation_large_grid(ii), 0);
    ratio_grid(ii) = f_grid(ii) / n ^ ( (beta(3)-1) / 2 );
    error_small_grid(ii) = abs(approximation_small_grid(ii) - f_grid(ii)) / f_grid(ii);
    error_large_grid(ii) = abs(approximation_large_grid(ii) - f_grid(ii)) / f_grid(ii);
end


% Main graph: where each approximation is good.
figure();
p = semilogx(n_grid, f_grid);
    p.Color = 'red';
    p.LineWidth = 2;
    p.LineStyle = '-';
hold on;
p = semilogx(n_grid, approximation_small_grid);
    p.Color = 'blue';
    p.LineWidth = 2;
    p.LineStyle = '--';
p = semilogx(n_grid, approximation_large_grid);
    p.Color = 'green';
    p.LineWidth = 2;
    p.LineStyle = ':';
p = semilogx(n_grid, n_grid'*0 + f_infinity);
    p.Color = 'black';
    p.LineWidth = 0.5;
    p.LineStyle = '-';

ylim([0, 7e-3]);
xlim([1e0, 1e11]);
    
set(gca, 'FontSize', 14);
pbaspect([1.618, 1, 1]);
legend("Production function", "Small n approximation", "Large n approximation", "Value of perfect information", ...
    "Location", "Northwest");
xlabel('Sample size $n$', 'interpreter', 'latex');
ylabel("Total product (% gain in success rate)");

print(gcf, '-depsc2', './output/figures/approximation-range.eps');