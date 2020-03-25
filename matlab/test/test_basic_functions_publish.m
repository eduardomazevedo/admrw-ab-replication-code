%% Start
addpath './classes';
addpath './functions';

%% Publish
publish_options = struct();
publish_options.format = 'html';
publish_options.outputDir = './test/test_basic_functions';
publish('./test/test_basic_functions.m', publish_options);
publish_options.format = 'pdf';
publish('./test/test_basic_functions.m', publish_options);