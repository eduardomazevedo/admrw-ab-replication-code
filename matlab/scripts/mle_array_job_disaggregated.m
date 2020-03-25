%% 1) Start
diary('./log/mle_array_job_disaggregated.m.log');
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
use_weights = 1; % Use weights or not.
algorithm = "gradient knitro";
knitro_flag = 1;


%% 3) Load data
data           = readtable('./data/intermediate/convolution-data-disaggregated.csv', "Delimiter", ",");


%% 4) Create variables
% Create a list of specifications to run through
n_specifications = 12;

if job_number > n_specifications
    error('Job number is greater than the number of specifications.');
end


%% 6) Estimation
% Start timer.
tic;


% Select subset of the data
rows = job_number == data.Spec;
rows = find(rows);
% rows = datasample(rows, 10); % Turn on for imprecise fast estimation.

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
[beta, l, flag, output, variance_matrix] = ...
    Twee.fit_g(...
    z, sigma, beta_initial, ...
    g, dg, dgg, ...
    max_iterations, ...
    knitro_flag, ...
    estimation_weights);

time = toc;


%% 7) Collect table of results
algorithm = {algorithm};
metric = {"session success rate"};
variance_matrix = {variance_matrix};
spec_number = job_number;
spec_description = data.SpecDescription(rows(1));
n_observations = length(legit_probability);
n_effective_observations = sum(legit_probability);

mle_output = table(spec_number, spec_description, ...
                   use_weights, algorithm, metric, ...
                   n_observations, n_effective_observations, ...
                   beta, variance_matrix, ...
                   l, ... 
                   time, ...
                   flag, output);
               

%% 5) End
filename = ['./matlab/mat/mle-disaggregated', '-', num2str(job_number), '.mat'];
save(filename, 'mle_output');
disp(['mle job #', num2str(job_number), ' done.']);
display(mle_output);
fprintf('\nAnalysis ended on %s\n\n', datetime);
diary off;