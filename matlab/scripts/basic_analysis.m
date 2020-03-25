%% basic analysis
% This script performs basic analyses and saves html output. It is for our
% internal use.


%% Start
close all
addpath './matlab/classes';
addpath './matlab/functions';


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

% Restrict attention to knitro estimate without weights.
mle_output = mle_output(strcmp(mle_output.algorithm, 'gradient knitro') & ...
    (mle_output.use_weights == 0), :);

row = strcmp(mle_output.metric, metric);
row = find(row);
beta = mle_output(row, 'beta');
beta = table2array(beta);


%% Set variables
data = data(strcmp(data.Metric, metric), :);
delta_hat = data.Delta;
sigma = data.StdErrorDelta;
sigma_n = data.sigma;

% Sort by z.
[delta_hat, sorted_index] = sort(delta_hat);
sigma = sigma(sorted_index);
sigma_n = sigma_n(sorted_index);


%% Print estimates
display(metric);
display(mle_output(:, {'metric', 'beta'}));


%% q-q plot of z
% TODO PEPE: are these axis labels switched?
qqrange = 0.01:0.01:0.99;
aux1 = quantile((delta_hat - mean(delta_hat)) ./ (std(delta_hat)), qqrange);
aux2 = norminv(qqrange, 0, 1);

figure();
    scatter(aux2, aux1);
    hold on;
    %45 degree line
        plot(aux2, aux2, 'r');
    hold off
    axis([min(aux1) max(aux1) min(aux1) max(aux2)]);
    xlabel('Quantiles of the Standard Normal');
    ylabel('Quantiles of the Standardized Data');


%% Compare the empirical c.d.f with the mixture c.d.f
figure()
    fig1 = cdfplot(delta_hat);
    set(fig1, 'LineWidth', 3)
    xlabel('z');
    ylabel('F(z)');
    hold on;
%   % TODO: This is not working right now. It is pretty easy but needs
%   fixing.
%     plot(z, ...
%         mean(Twee.L_cdf(z, sigma, ...
%         beta, g), 2), '--r');
    hold off;
    legend('Empirical Distribution of z','Estimated Distribution of z')
    legend('location','southeast')
    legend boxoff


%% Plot tweedie estimate of posterior mean vs signal.
Z = -0.25:0.01:0.25;
figure();
    plot(Z, Twee.mean_posterior(Z, mean(sigma), beta, g));
    title('Tweedie estimate vs signal');
    
figure();
    precision_noise  = mean(data.StdErrorDelta)^(-2);
    precision_signal = beta(2)^(-2);
    normal_posterior_shrinkage = precision_signal ./ (precision_signal + precision_noise);
    normal_posterior = normal_posterior_shrinkage .* Z + ...
        (1-normal_posterior_shrinkage) .* beta(1);
    plot(Z, Twee.mean_posterior(Z, mean(sigma), beta, g));
    hold on;
    plot(Z, normal_posterior);
    title('Tweedie estimate and Gaussian estimate vs. signal');
    
    
%% Business conversion from ssr to yearly revenue
    tousdbn = 10 / 100 * 6;
    tousdM  = tousdbn * 1e3;
    tousd   = tousdbn * 1e9;
    parallel_experiments_per_user = 1; % I will just make this equal to one for now, and it is easy to convert if we want to make an assumption about it.
    weeks_per_year = 52.14;


%% Decomposition of gains.
% Define useful variables
posterior = Twee.mean_posterior(delta_hat, sigma, beta, g);
a_p_value = delta_hat ./ sigma > 1.96;
a_z_positive = delta_hat > 0;
a_bayesian = posterior > 0;
percentile_vector = (1:length(delta_hat)) ./ length(delta_hat);
probability_legit = data.ProbabilityLegitFitted(sorted_index);

% 1. p value policy with incorrect frequentist measure and Bayesian measure.
figure();
    plot(percentile_vector, [...
        cumsum(delta_hat .* a_p_value), ...
        cumsum(posterior .* a_p_value)]);
    xlim([0, 1.01])
    legend({'Frequentist Estimate', 'Bayesian Estimate'});
    title('Decomposition of past gains: status quo policy');

if strcmp(metric, 'session success rate')
    fprintf('\n\nEffectiveness of the p-value policy:\n');
    fprintf('Naive effectiveness from %d ideas, out of which %d are implemented, naively looks like a gain of %.2f%% = %.0f million USD per year.\n', ...
        round(sum(probability_legit)), ...
        round(sum(a_p_value .* probability_legit)), ...
        sum(delta_hat .* a_p_value .* probability_legit), ...
        sum(delta_hat .* a_p_value .* probability_legit) * tousdM);
    fprintf('Bayesian gain is %.2f%% = %.0f million USD per year.\n', ...
        sum(posterior .* a_p_value .* probability_legit), ...
        sum(posterior .* a_p_value .* probability_legit) * tousdM);
    fprintf('The marginal idea with a p value of 5%% has a gain of %.3f%% = %.2f million USD per year.\n', ...
        min(posterior(a_p_value)), min(posterior(a_p_value)) * tousdM);
end

% 2. z > 0 policy with incorrect frequentist measure and Bayesian measure.
figure();
    plot(percentile_vector, [...
        cumsum(delta_hat .* a_z_positive .* probability_legit), ...
        cumsum(posterior .* a_z_positive .* probability_legit)]);
    xlim([0, 1.01])
    legend({'Frequentist Estimate', 'Bayesian Estimate'});
    title('Decomposition of past gains: $z>0$ policy');

% 3. The three policies evaluated correctly by the Bayesian method.
figure();
    plot(percentile_vector, [...
        cumsum(posterior .* a_p_value .* probability_legit), ...
        cumsum(posterior .* a_z_positive .* probability_legit), ...
        cumsum(posterior .* a_bayesian .* probability_legit)]);
    xlim([0, 1.01])
    legend({'p-value', 'z>0', 'Bayesian'});
    title('Real past gains by policy');

if strcmp(metric, 'session success rate')
    fprintf('\n\nEffectiveness of Bayesian policy: implementing more ideas with small gains:\n');
    fprintf('Naive effectiveness from %d ideas, out of which %d are implemented, naively looks like a gain of %.2f%% = %.0f million USD per year.\n', ...
        round(sum(probability_legit)), ...
        round(sum(a_bayesian .* probability_legit)), ...
        sum(delta_hat .* a_bayesian .* probability_legit), ...
        sum(delta_hat .* a_bayesian .* probability_legit) * tousdM);
    fprintf('Bayesian gain is %.2f%% = %.0f million USD per year.\n', ...
        sum(posterior .* a_bayesian .* probability_legit), ...
        sum(posterior .* a_bayesian .* probability_legit) * tousdM);
    fprintf('We gain an extra %.2f%% = %.0f million USD per year at the cost of implementing an extra %d ideas.\n', ...
        sum(posterior .* (a_bayesian - a_p_value) .* probability_legit), ...
        sum(posterior .* (a_bayesian - a_p_value) .* probability_legit) * tousdM, ...
        round(sum((a_bayesian - a_p_value) .* probability_legit)));  
    fprintf('We can use the Bayesian estimator to figure out what is a good tradeoff, as opposed to relying on p value.\n');
    fprintf('If we require an improvement of 0.003%% = 1.8M USD per year, the gain is %.0f million USD from implementing %d ideas.\n', ...
        sum(posterior .* (posterior >= 0.003) .* probability_legit) * tousdM, ...
        round(sum((posterior >= 0.003) .* probability_legit)));
    fprintf('Relative to the p-value policy, this is an extra %.0f million USD per year at the cost of implementing %d ideas.\n', ...
        sum(posterior .* ((posterior >= 0.003) - a_p_value) .* probability_legit) * tousdM, ...
        round(sum(((posterior >= 0.003) - a_p_value) .* probability_legit)));
end

% 4. Fraction of the gains coming from outliers
if strcmp(metric, 'session success rate')
    t_stats = delta_hat./sigma;
    total_gain_frequentist = sum(posterior .* a_p_value .* probability_legit);
    fprintf('\n\nImportance of outliers in the p-value policy:\n');
    fprintf('%.2f%% of the gains comes from innovations 3*se\n', ...
        100*sum(posterior .* (t_stats > 3) .* probability_legit)/total_gain_frequentist);
    fprintf('%.2f%% of the gains comes from innovations 4*se\n', ...
        100*sum(posterior .* (t_stats > 4) .* probability_legit)/total_gain_frequentist);
end

% 5. Fraction of the gain coming from the top few obs
if strcmp(metric, 'session success rate')
    q = 0.05;
    gain = posterior .* a_p_value .* (delta_hat >= quantile(delta_hat, 1 - q)) .* probability_legit;
    gain = sum(gain);
    fprintf('%.2f%% of the gains of the p value policy come from the top 5%% innovations\n', ...
        100 * gain / total_gain_frequentist);
    q = 0.02;
    gain = posterior .* a_p_value .* (delta_hat >= quantile(delta_hat, 1 - q)) .* probability_legit;
    gain = sum(gain);
    fprintf('%.2f%% of the gains of the p value policy come from the top 2%% innovations\n', ...
        100 * gain / total_gain_frequentist);
    q = 0.01;
    gain = posterior .* a_p_value .* (delta_hat >= quantile(delta_hat, 1 - q)) .* probability_legit;
    gain = sum(gain);
    fprintf('%.2f%% of the gains of the p value policy come from the top 1%% innovations\n', ...
        100 * gain / total_gain_frequentist);
end



%% Production function
%Set up a value of sigma_i
mean_sigma_individual = mean(sigma_n);
n_grid = linspace(0, 30*10^6, n_points_production_function);
Y = zeros(3, n_points_production_function);

%Estimate the production function using the Monte-Carlo draws of delta and
%epsilon

for i = 1:n_points_production_function
    % optimal policy
    Y(1,i) = ...
        Twee.f(n_grid(i), mean_sigma_individual, beta, g);
    % z > 0 policy
    Y(2,i) = ...
        Twee.f(n_grid(i), mean_sigma_individual, beta, g, 0);
    % p value policy
    z_bar = 1.96 .* mean_sigma_individual ./ sqrt(n_grid(i));
    Y(3,i) = ...
        Twee.f(n_grid(i), mean_sigma_individual, beta, g, z_bar);
end

% Production function
figure();
    plot(n_grid', Y');
    grid on;
    title('Production Function under different policies');
    xlabel('Size of the Experiment (n)');
    legend(...
        'optimal',...
        'z > 0',...
        'p value');
    legend('Location','Southeast');
    
if strcmp(metric, 'session success rate')
    fprintf('\n\nUsing the production function to value data.\n');
    fprintf('\n\nUsing the production function to value data.\n');    
end

% Production function: log x axis.
figure();
    plot(log10(n_grid'), Y');
    grid on;
    title('Production Function under different policies -- loged x axis');
    xlabel('log Size of the Experiment (n)');
    legend(...
        'optimal',...
        'z > 0',...
        'p value');
    legend('Location','Southeast');
    legend boxoff;

% Average and marginal product
AP = Y(1, :) ./ n_grid;
MP = gradient(Y(1,:), n_grid);
figure();
    plot(n_grid', [AP; MP]');
    grid on;
    title('Average and Marginal Product');
    xlabel('Size of the Experiment (n)');
    legend(...
        'AP',...
        'MP');
    legend('Location','Southeast');
    legend boxoff;

% Average and marginal product: log-log plot
figure();
    plot(log10(n_grid'), log10([AP; MP]'));
    grid on;
    title('Average and Marginal Product -- log-log plot');
    xlabel('log Size of the Experiment (n)');
    legend(...
        'AP',...
        'MP');
    legend('Location','Southeast');
    legend boxoff;

% Average and marginal product ratio
figure();
    plot(n_grid', MP' ./ AP');
    grid on;
    title('Marginal to average product ratio');
    xlabel('Size of the Experiment (n)');
    legend('Location','Southeast');
    legend boxoff;

% Robustness of the production function with respect to the tail
% coefficient
YY = zeros(3, n_points_production_function);
for i = 1:n_points_production_function
    YY(1, i) = Twee.f(n_grid(i), mean_sigma_individual, beta + [0 0 1], g);
    YY(2, i) = Twee.f(n_grid(i), mean_sigma_individual, beta + [0 0 2], g);
    YY(3, i) = Twee.f(n_grid(i), mean_sigma_individual, beta + [0 0 3], g);
end
figure();
    plot(n_grid', [Y(1, :); YY]');
    grid on;
    title('Production Function with different tail coefficients');
    xlabel('Size of the Experiment (n)');
    legend(...
        'estimated',...
        'tail + 1',...
        'tail + 2',...
        'tail + 3');
    legend('Location','Southeast');
    
    
%% Production function interpretation and dollar values
if strcmp(metric, 'session success rate')
    ap20 = Twee.f(20e6, mean_sigma_individual, beta, g) / 20e6;
    ap5  = Twee.f(5e6 , mean_sigma_individual, beta, g) /  5e6;
    ap1  = Twee.f(1e6 , mean_sigma_individual, beta, g) /  1e6;
    
    mp20 = interp1q(n_grid', MP', 20e6);

    % User-weeks
    fprintf('\n\nValue of data\n');
    fprintf('Assuming that 1%% SSR is worth 10%% yearly revenue, and yearly revenue of 6bn (following EXP slides)\n');
    fprintf('Assuming %d parallel experiments\n', parallel_experiments_per_user);
    fprintf('The gains are in yearly USD as opposed to once off.\n');
    fprintf('Average product of data is %.2f dollars/1,000 user-week in 20M users experiment, %.2f in 5M user experiment, and %.2f in 1M user experiment.\n', ...
        ap20 * tousd * parallel_experiments_per_user * 1e3, ...
        ap5 * tousd * parallel_experiments_per_user * 1e3, ...
        ap1 * tousd * parallel_experiments_per_user * 1e3);
    fprintf('Marginal value of data in 20M user experiment is only %.2f dollars/1,000 user-week\n', mp20 * tousd * parallel_experiments_per_user * 1e3);

    % User-years
    fprintf('\n\nValue of data\n');
    fprintf('Assuming that 1%% SSR is worth 10%% yearly revenue, and yearly revenue of 6bn (following EXP slides)\n');
    fprintf('Assuming %d parallel experiments\n', parallel_experiments_per_user);
    fprintf('The gains are in yearly USD as opposed to once off.\n');
    fprintf('Average product of data is %.2f dollars/user-year in 20M users experiment, %.2f in 5M user experiment, and %.2f in 1M user experiment.\n', ...
        ap20 * tousd * parallel_experiments_per_user * weeks_per_year, ...
        ap5 * tousd * parallel_experiments_per_user * weeks_per_year, ...
        ap1 * tousd * parallel_experiments_per_user * weeks_per_year);
    fprintf('Marginal value of data in 20M user experiment is only %.2f dollars/ user-year\n', ...
        mp20 * tousd * parallel_experiments_per_user * weeks_per_year);    
    
    
    speed_up_10pc = 1.1 * Twee.f(20e6 / 1.1, mean_sigma_individual, beta, g) / ...
        Twee.f(20e6, mean_sigma_individual, beta, g);
    speed_up_20pc = 1.2 * Twee.f(20e6 / 1.2, mean_sigma_individual, beta, g) / ...
        Twee.f(20e6, mean_sigma_individual, beta, g);    
    speed_up_doubling = 2 * Twee.f(20e6 / 2, mean_sigma_individual, beta, g) / ...
        Twee.f(20e6, mean_sigma_individual, beta, g);  
    
    fprintf('\n');
    fprintf('Increasing the number of experiments by 10%% speeds innovation by %.2f%%\n', 100*(speed_up_10pc-1));
    fprintf('Increasing the number of experiments by 20%% speeds innovation by %.2f%%\n', 100*(speed_up_20pc-1));
    fprintf('Doubling the number of experiments speeds innovation by %.2f%%\n', 100*(speed_up_doubling-1));
end

    
%% Production function with small n
%Set up a value of sigma_i
n_points_production_function = 40;
mean_sigma_individual = mean(sigma_n);
n_grid_small = linspace(0, 10^4, n_points_production_function);
Y = zeros(4, n_points_production_function);

%Estimate the production function using the Monte-Carlo draws of delta and
%epsilon

for i = 1:n_points_production_function
    % regular
    Y(1, i) = ...
        Twee.f(n_grid_small(i), mean_sigma_individual, beta + [0 0 0], g);
    % tail + 1
    Y(2, i) = ...
        Twee.f(n_grid_small(i), mean_sigma_individual, beta + [0 0 1], g);
    % tail + 2
    Y(3, i) = ...
        Twee.f(n_grid_small(i), mean_sigma_individual, beta + [0 0 2], g);
    % tail + 3
    Y(4, i) = ...
        Twee.f(n_grid_small(i), mean_sigma_individual, beta + [0 0 3], g);
end

figure();
    plot(n_grid_small', Y');
    grid on;
    title('Production Function with different tail coefficients - small n');
    xlabel('Size of the Experiment (n)');
    legend(...
        'estimated',...
        'tail + 1',...
        'tail + 2',...
        'tail + 3');
    legend('Location','Southeast');


%% End
fprintf('\nAnalysis ended on %s\n\n', datetime);
diary off;
