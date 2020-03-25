%% Generates a table of disaggregated mle results.
clear;
addpath './matlab/functions';
load('./matlab/mat/mle_disaggregated.mat');
mle_output_disaggregated = sortrows(mle_output_disaggregated, 'spec_number');

% Create a matrix with the data for the table
% Format:
% specification, M, s, alpha, observations headers.
% Values in odd lines. Stars for alpha < 3: * = 99%. ** = 99.9%
% Standard errors in even lines right after
tostringsci = @(ss) num2str(ss, '%0.2e');
tostringreg = @(ss) num2str(ss, '%0.2f');
n_specs = height(mle_output_disaggregated);

tab = table();
for ii = 1:n_specs
    M = mle_output_disaggregated.beta(ii, 1);
    s = mle_output_disaggregated.beta(ii, 2);
    alpha = mle_output_disaggregated.beta(ii, 3);
    n_observations = mle_output_disaggregated.n_observations(ii);
    covariance_matrix = mle_output_disaggregated.variance_matrix{ii};
    se_vector = sqrt(diag(covariance_matrix));
    t_stat_alpha = (3 - alpha) / se_vector(3);
    if t_stat_alpha > norminv(0.999)
        stars = '***';
    elseif t_stat_alpha > norminv(0.99)
        stars = '**';
    elseif t_stat_alpha > norminv(0.95)
        stars = '*';    
    else
        stars = '';
    end

    % row with M, s, alpha, observations
    row = table( ...
        mle_output_disaggregated.spec_description(ii), ...
        {tostringsci(M)}, {tostringsci(s)}, {[tostringreg(alpha), stars]}, {num2str(n_observations)});
    tab = [tab; row];

    % row with blank, standard errors, blank
    row = table( ...
        {''}, ...
        {['(', tostringsci(se_vector(1)), ')']}, {['(', tostringsci(se_vector(2)), ')']}, {['(', tostringreg(se_vector(3)), ')']}, {''});
    tab = [tab; row];

end

% Generate the table using the latexTable function
input = struct();
input.data = tab;
input.tableBorders              = 0;
input.tableColLabels = {'Subsample', '$M$','$s$','$\alpha$', 'observations'};

output = latexTable(input);
output = output(5:length(output) - 4); % Hack to only get the tabular environment, not table. We also want to make our own header.

% Save output
fid = fopen("./output/tables/mle-disaggregated-2.tex", "w");
    fprintf(fid,'%s\n','\begin{tabular}{lcccc}');
    fprintf(fid,'%s\n','\toprule');
    fprintf(fid,'%s\n',strcat('Subsample', '&', '$M$', '&', '$s$', '&', '$\alpha$', '&', 'observations', '\\'));
    fprintf(fid,'%s\n','\midrule');
    fprintf(fid, '%s\n', output{:});
    fprintf(fid,'%s\n','\bottomrule');
    fprintf(fid,'%s\n','\end{tabular}');
fclose(fid);
