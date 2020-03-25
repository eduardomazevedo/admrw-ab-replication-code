%% 1) Start
diary('./log/mle_array_job.m.log');
fprintf('\nAnalysis started on %s\n\n', datetime);
clear;
close all;
addpath './matlab/classes';
addpath './matlab/functions';
job_number = str2double(getenv('SGE_TASK_ID'));


%% 2) Parameters
g              = @t_distribution;    % Prior
dg             = @d_t_distribution;  % First derivative of the prior
dgg            = @dd_t_distribution; % Second derivative of the prior
max_iterations = 200;                % Iteration limits

% We will do several specifications varying whether we use weights and
% varying the algorithm.
specification_use_weights = {0, 1}'; % Use weights or not.
specification_algorithm = {'no gradient', 'gradient fmincon', 'gradient knitro'}';


%% 3) Load data
data           = readtable('./data/intermediate/convolution-data.csv');
specification_metrics        = unique(data.Metric);


%% 4) Create variables
% Create a list of specifications to run through
n_specifications = ...
    length(specification_use_weights) * ...
    length(specification_algorithm) * ...
    length(specification_metrics);

if job_number > n_specifications
    error('Job number is greater than the number of specifications.');
end

specification_list = cell(n_specifications, 1);

ii = 1;
for i1 = 1:length(specification_use_weights)
    for i2 = 1:length(specification_algorithm)
        for i3 = 1:length(specification_metrics)
            specification_list{ii} = ...
                {...
                specification_use_weights{i1}, ...
                specification_algorithm{i2}, ...
                specification_metrics{i3}, ...
                };
            ii = ii + 1;
        end
    end
end


%% 6) Estimation
% Start timer.
tic;

% Set specification parameters
specification = specification_list{job_number};
use_weights = specification{1};
algorithm = specification{2};
metric = specification{3};

if strcmp(algorithm, 'gradient knitro')
    knitro_flag = 1;
else
    knitro_flag = 0;
end

% Select subset of the data
rows = strcmp(metric, data.Metric);
rows = find(rows);
%rows = datasample(rows, 10); % Turn on for imprecise fast estimation.

% Create estimation data
z     = data.Delta(rows);
sigma = data.StdErrorDelta(rows);    % bad notation.
legit_probability = data.ProbabilityLegitFitted(rows);

% Initial value of beta
beta_initial      = [mean(z), std(z), 3];

% Weights
if use_weights
    estimation_weights = legit_probability;
else
    estimation_weights = ones(length(rows), 1);
end

% Estimation
if strcmp(algorithm, 'gradient fmincon') || ...
        strcmp(algorithm, 'gradient knitro')
    [beta, l, flag, output, variance_matrix] = ...
        Twee.fit_g(...
        z, sigma, beta_initial, ...
        g, dg, dgg, ...
        max_iterations, ...
        knitro_flag, ...
        estimation_weights);

elseif strcmp(algorithm, 'no gradient')
    % TODO: the no gradient method does not allow for weights!
    [beta, l, flag, output] = ...
        Twee.fit_g_no_gradient(...
        z, sigma, beta_initial, ...
        g, max_iterations);
    variance_matrix = nan(3, 3);
end

time = toc;


%% 7) Collect table of results
algorithm = {algorithm};
metric = {metric};
variance_matrix = {variance_matrix};
mle_output = table(use_weights, algorithm, metric, ...
                   beta, variance_matrix, ...
                   l, ... 
                   time, ...
                   flag, output);
               

%% 5) End
filename = ['./matlab/mat/mle', '-', num2str(job_number), '.mat'];
save(filename, 'mle_output');
disp(['mle job #', num2str(job_number), ' done.']);
display(mle_output);
fprintf('\nAnalysis ended on %s\n\n', datetime);
diary off;