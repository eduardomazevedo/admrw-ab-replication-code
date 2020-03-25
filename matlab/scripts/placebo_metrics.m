% Additional results requested by JPE reviewers in reject and resubmit
% Plot the production function and posterior mean function for all metrics.
% The goal is to show that we recover the common sense result that the long
% run metrics are not actionable, while the short run metrics are all
% qualitatively similar.
%% Start
% 
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

% Select data for success rate
% data = data(strcmp(data.Metric, 'session success rate'), :);

% Select MLE output with preferred specification of no weights.
row = (mle_output.use_weights == 0);
mle_output = mle_output(row, :);

% Choose beta for success rate with benchmark estimate
% row = strcmp(mle_output.metric, 'session success rate');
% beta = mle_output(row, 'beta');
% beta = table2array(beta);

% Create useful variables
% delta_hat = data.Delta;
% sigma_n = data.StdErrorDelta;
% sigma = data.sigma;

% mean_sigma = mean(sigma);

% Sort variables by the signal
% [delta_hat, sorted_index] = sort(delta_hat);
% sigma = sigma(sorted_index);
% sigma_n = sigma_n(sorted_index);


%% Fix metric names for paper
mle_output.metric = strrep(mle_output.metric, 'revenue', 'na');
mle_output.metric = strrep(mle_output.metric, 'utility', 'na');
mle_output.metric = strrep(mle_output.metric, 'page load time', 'na');
mle_output.metric = strrep(mle_output.metric, 'session success rate', 'session success rate');
mle_output.metric = strrep(mle_output.metric, 'page click rate', 'short-run metric 1');
mle_output.metric = strrep(mle_output.metric, 'quickback rate', 'na');
mle_output.metric = strrep(mle_output.metric, 'time to success', 'short-run metric 3');
mle_output.metric = strrep(mle_output.metric, 'queries per user', 'long-run metric 1');
mle_output.metric = strrep(mle_output.metric, 'sessions per user', 'long-run metric 2');

row = ~strcmp(mle_output.metric, 'na');
mle_output = mle_output(row, :);
mle_output = sortrows(mle_output, 'metric');

data.Metric = strrep(data.Metric, 'revenue', 'na');
data.Metric = strrep(data.Metric, 'utility', 'na');
data.Metric = strrep(data.Metric, 'page load time', 'na');
data.Metric = strrep(data.Metric, 'session success rate', 'session success rate');
data.Metric = strrep(data.Metric, 'page click rate', 'short-run metric 1');
data.Metric = strrep(data.Metric, 'quickback rate', 'na');
data.Metric = strrep(data.Metric, 'time to success', 'short-run metric 3');
data.Metric = strrep(data.Metric, 'queries per user', 'long-run metric 1');
data.Metric = strrep(data.Metric, 'sessions per user', 'long-run metric 2');

row = ~strcmp(data.Metric, 'na');
data = data(row, :);
data = sortrows(data, 'Metric');

metrics_list = mle_output.metric;
n_metrics = length(metrics_list);


%% Make calculations: production function and posterior mean
% Define grid
n_points = 100;
y = zeros(1, n_metrics);
p_hat = zeros(n_points, n_metrics);
f_infinity = zeros(1, n_metrics);
delta_grid = linspace(-0.6, 0.6, n_points);

for ii = 1:n_metrics
    metric = metrics_list(ii);
    
    row_mle = strcmp(metric, mle_output.metric);
    beta = mle_output(row_mle, 'beta');
    beta = table2array(beta);
    
    row_data = strcmp(metric, data.Metric);
    sigma = mean(table2array(data(row_data, 'sigma')));
    largest_delta_hat = max(abs(table2array(data(row_data, 'Delta'))));
    interquartile_range = iqr(table2array(data(row_data, 'Delta')));
    
    f_infinity(ii) = Twee.f(Inf, sigma, beta, g);
    
    y(ii) = Twee.f(20e6, sigma, beta, g) / f_infinity(ii);
    
    for jj = 1:n_points
        delta = delta_grid(jj);
        p_hat(jj, ii) = Twee.mean_posterior(delta, sigma ./ sqrt(20e6), beta, g);
    end
end


%% Plot graphs
% p hat
figure();
set(gca, 'FontSize', 18);
hold on;

p = plot(delta_grid, p_hat(:, 1));
    p.Color = 'red';
    p.LineWidth = 0.5;
    p.LineStyle = '-';
p = plot(delta_grid, p_hat(:, 2));
    p.Color = 'red';
    p.LineWidth = 0.5;
    p.LineStyle = '--';
p = plot(delta_grid, p_hat(:, 3));
    p.Color = 'black';
    p.LineWidth = 1.5;
    p.LineStyle = '-';
p = plot(delta_grid, p_hat(:, 4));
    p.Color = 'black';
    p.LineWidth = 1.5;
    p.LineStyle = '--';
p = plot(delta_grid, p_hat(:, 5));
    p.Color = 'black';
    p.LineWidth = 1.5;
    p.LineStyle = ':';
    
legend(metrics_list, "Location", "Northwest")

% Labels
xlabel({'Signal $\widehat{\delta} _i$'},'Interpreter','latex');
ylabel({'Posterior mean $P_i(\widehat{\delta} _i, n_i)$'},'Interpreter','latex');

% Save
print(gcf,'-depsc2','./output/figures/placebos-posterior-mean.eps');


% Production function
figure();
set(gca, 'FontSize', 18);

bar(categorical(metrics_list), 100 * y);
ylim([0, 100]);
ytickformat('percentage');
for ii = 1:2
    text(ii, 100 * y(ii), ...
        [num2str(y(ii), '%0.2e'), '%'], ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'bottom')
end
for ii = 3:n_metrics
    text(ii, 100 * y(ii), ...
        [num2str(100 * y(ii), '%0.1f'), '%'], ...
               'HorizontalAlignment', 'center', ...
               'VerticalAlignment', 'bottom')
end

ylabel('Percentage of gains obtained in 20 million user experiment $f_i(20e6)/f_i(\infty)$', 'interpreter', 'latex');
% Save
print(gcf,'-depsc2','./output/figures/placebos-production.eps');
