%% Start
clear;
close all
addpath './matlab/classes';
addpath './matlab/functions';
metric_i = 1;


%% Parameters
% Prior
g  = @t_distribution;
dg = @d_t_distribution;

% Computation
n_points_production_function = 400;


%% Load data
data = readtable('./data/intermediate/convolution-data.csv');
metrics = unique(data.Metric);
metric = metrics{metric_i};
metric_no_space = strrep(metric, ' ', '_');

diary(['./log/basic-analysis/' metric_no_space '.log']);
fprintf('\nAnalysis started on %s\n\n', datetime);

mle_output = load('./matlab/mat/mle.mat');
mle_output = mle_output.mle_output;
beta = mle_output(metric_i, 'beta');
beta = table2array(beta);

% Filter
% data = data(data.NTreatments == 1 & data.AnalysisRange == 7, :);

% Calculate summary stats
aux = grpstats(data, 'Metric', 'mean', 'DataVars', 'StdErrorDelta');
mle_output.mean_se = aux.mean_StdErrorDelta;
mle_output.theta = mle_output.beta(:, 2).^2 ./ (mle_output.beta(:, 2).^2 + mle_output.mean_se.^2);


%% Set variables
data = data(strcmp(data.Metric, metric), :);
z = data.Delta;
sigma = data.StdErrorDelta;
sigma_individual = data.sigma;

% Sort by z.
[z, sorted_index] = sort(z);
sigma = sigma(sorted_index);
sigma_individual = sigma_individual(sorted_index);

% Set probability legit
probability_legit = data.ProbabilityLegitFitted(sorted_index);


%% Print estimates
display(metric);
display(mle_output(:, {'metrics', 'beta', 'mean_se', 'theta'}));


%% Variables for graph
% Basic
mean_sigma_individual = mean(sigma_individual);
mean_sigma_experiment = mean_sigma_individual / sqrt(20e6);

% Grid
leftmost_z = -0.05;
z_grid = linspace(leftmost_z, 0.4, 100);

% Posterior
% mean_posterior_sample = Twee.mean_posterior(z, mean_sigma_experiment, beta, g);
N = length(z);
n = floor(0.98 .* N);
z_tail = z(n : N);
z_head = z(1 : n - 1);
mean_posterior_sample_tail = mean_posterior_sample(n : N);
mean_posterior_sample_head = mean_posterior_sample(1 : n-1);

% Density
L = @(zz) Twee.L(zz, mean_sigma_experiment, beta, g);
L_grid = L(z_grid);

%% Make graph
close all;
fig = figure();

set(fig, ...
    'defaultAxesColorOrder', [0 0 0; 232/255 76/255 79/255]);
pbaspect([16 9 9]);
xlabel('Percentage Improvement');

% Density
yyaxis left;
plot(z_grid, L_grid, ...
    'Color', [0.5 0.5 0.5]);

ax.YColor = 'black';

% Dots and horizontal line
marker_size = 35;
yyaxis right;
hold on;
ax.YColor = [232,76,79]/255;
ylabel('Percentage Improvement');

hline = refline([0 0]);
hline.Color = [232,76,79]/255;

scatter(z_head, mean_posterior_sample_head, ...
    'MarkerFaceColor', 'black', ...
    'MarkerEdgeColor', 'black', ...
    'MarkerFaceAlpha', .017, ...
    'MarkerEdgeAlpha', .017, ...
    'SizeData', marker_size);
% Scatter version
scatter(z_tail, mean_posterior_sample_tail, ...
    'MarkerFaceColor', [232,76,79]/255, ...
    'MarkerEdgeColor', [232,76,79]/255, ...
    'MarkerFaceAlpha', 1, ...
    'MarkerEdgeAlpha', 1, ...
    'SizeData', marker_size);
% Stem version
% stem(z_tail, mean_posterior_sample_tail, ...
%     'MarkerFaceColor', [232,76,79]/255, ...
%     'MarkerEdgeColor', [232,76,79]/255, ...
%     'LineStyle', '-');

% Labels
% text(0.121892291, 0.121892291 - 0.025, ...
%     'New Header', ...
%     'Color' , [232,76,79]/255);
% text(0.13087864 + 0.025, 0.13087864, ...
%     'Local Search Improvement', ...
%     'Color' , [232,76,79]/255);

axis([leftmost_z 0.4 leftmost_z, 0.4]);