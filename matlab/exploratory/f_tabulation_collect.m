%% 1) Start
clear;
close all;
addpath './matlab/classes';
addpath './matlab/functions';


%% 2) Collect output data
files = dir('./matlab/mat/tabulate-f-*.mat');
files = {files.name};

% Create variables
load(['./matlab/mat/' files{1}]);
y_output = nan(length(grid_n), length(grid_alpha));

for file_name = files
    load(['./matlab/mat/' file_name{1}]);
    y_output(job_number, :) = y;
end

y = y_output;
clear y_output files job_number n_files;

save('./matlab/mat/tabulate-f.mat')