%% Test MLE with simulated data
%% Start
addpath ./matlab/classes;
addpath ./matlab/functions;
clear;
close all;

% Print date
fprintf('\nAnalysis started on %s\n\n', datetime);


%% Parameters
% Simulation
n_draws = 1e2;

% Estimation
g              = @t_distribution;    % Prior
dg             = @d_t_distribution;  % First derivative of the prior
dgg            = @dd_t_distribution; % Second derivative of the prior

max_iterations = 200;                % Iteration limits
knitro_flag = 0;

sigma_n = 1;


%% Noise terms
epsilon = sigma_n * randn(n_draws, 1);


%% Normal distribution estimated as a t distribution
% Simulation
delta = randn(n_draws, 1);
mle_weights = ones(n_draws, 1);

delta_hat = delta + epsilon;

% Estimation
beta_initial = [mean(delta_hat), std(delta_hat), 3];

[beta_normal, l, flag, output, variance_matrix] = ...
        Twee.fit_g_conf(...
        delta_hat, sigma_n, beta_initial, ...
        g, dg, dgg, ...
        max_iterations, ...
        knitro_flag, mle_weights);
    
    
%% t distribution estimated as a t distribution
% Simulation
alpha = 1.5;
delta = trnd(alpha, [n_draws, 1]);
mle_weights = ones(n_draws, 1);

delta_hat = delta + epsilon;

% Estimation
beta_initial = [mean(delta_hat), std(delta_hat), 3];

[beta_t, l, flag, output, variance_matrix] = ...
        Twee.fit_g_conf(...
        delta_hat, sigma_n, beta_initial, ...
        g, dg, dgg, ...
        max_iterations, ...
        knitro_flag, mle_weights);