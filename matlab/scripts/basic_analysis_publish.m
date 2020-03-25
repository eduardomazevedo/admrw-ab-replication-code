%% Start
addpath './matlab/classes';
addpath './matlab/functions';
metric_i = str2double(getenv('SGE_TASK_ID'));


%% Get metric names
data = readtable('./data/intermediate/convolution-data.csv');
metrics = unique(data.Metric);
metric = metrics{metric_i};
metric_no_space = strrep(metric, ' ', '_');
clear data metrics metric;


%% Publish
publish_options = struct();
publish_options.format = 'html';
publish_options.outputDir = ['./output/basic_analysis/' metric_no_space];
publish('./matlab/scripts/basic_analysis.m', publish_options);
publish_options.format = 'pdf';
publish('./matlab/scripts/basic_analysis.m', publish_options);