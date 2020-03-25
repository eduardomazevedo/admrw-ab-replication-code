%% Start
diary('./log/f_tabulation_array_job.m.log');
fprintf('\nAnalysis started on %s\n\n', datetime);
clear;
close all;
addpath './matlab/classes';
addpath './matlab/functions';
job_number = str2double(getenv('SGE_TASK_ID'));


%% Parameters
n_points_alpha = 50;
n_points_n = 400;

left_n = 0;
right_n = 40e6;
left_alpha = 1.2;
right_alpha = 5;

g              = @t_distribution;    % Prior
dg             = @d_t_distribution;  % First derivative of the prior
dgg            = @dd_t_distribution; % Second derivative of the prior
max_iterations = 200;                % Iteration limits

sigma = 100;


%% Load data
mle_output = load('./matlab/mat/mle.mat');
mle_output = mle_output.mle_output;
row = ...
    (mle_output.use_weights == 0) & ...
    strcmp(mle_output.metric, 'session success rate') & ...
    strcmp(mle_output.algorithm, 'gradient knitro');
mle_output = mle_output(row, :);


%% Set variables
beta_estimate = mle_output.beta;

grid_n = linspace(left_n, right_n, n_points_n);
grid_alpha = linspace(left_alpha, right_alpha, n_points_alpha);

y = zeros(1, n_points_alpha);


%% Estimate the production function
for ii = 1:n_points_alpha
    % Parameters
    beta = beta_estimate;
    beta(3) = grid_alpha(ii);
    n = grid_n(job_number);
    
    y(ii) = Twee.f(n, sigma, beta, g);
end


%% Save
file_name = ['./matlab/mat/tabulate-f-' num2str(job_number) '.mat'];
save(file_name, 'y', 'grid_n', 'grid_alpha', 'job_number');