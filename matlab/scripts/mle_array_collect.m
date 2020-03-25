%% 1) Start
diary('./log/mle_array_collect.m.log');
fprintf('\nAnalysis started on %s\n\n', datetime);
clear;
close all;
addpath './matlab/classes';
addpath './matlab/functions';


%% 2) Collect output data
% Get file names
files = dir('./matlab/mat/mle-*.mat');
files = {files.name};
files_disaggregated = dir('./matlab/mat/mle-disaggregated-*.mat');
files_disaggregated = {files_disaggregated.name};
files = setdiff(files, files_disaggregated);
n = length(files);
n_disaggregated = length(files_disaggregated);

% Collect baseline output
rows = cell(n, 1);
for ii = 1:n
    rows{ii} = load(['./matlab/mat/', files{ii}]);
    rows{ii} = rows{ii}.mle_output;
    rows{ii}.output = {rows{ii}.output}; %Output has to be a cell array to be concatenated into the final table.
end

mle_output = vertcat(rows{:});

% Collect disaggregated output
rows = cell(n_disaggregated, 1);
for ii = 1:n_disaggregated
    rows{ii} = load(['./matlab/mat/', files_disaggregated{ii}]);
    rows{ii} = rows{ii}.mle_output;
    rows{ii}.output = {rows{ii}.output}; %Output has to be a cell array to be concatenated into the final table.
end

mle_output_disaggregated = vertcat(rows{:});

%% 5) End
save('./matlab/mat/mle.mat', 'mle_output');
save('./matlab/mat/mle_disaggregated.mat', 'mle_output_disaggregated');
disp('mle_array_collect.m done.');
display(mle_output);
display(mle_output_disaggregated);
fprintf('\nAnalysis ended on %s\n\n', datetime);
diary off;