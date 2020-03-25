%% This script generates the figures and constants for the paper


%% Start
diary('./log/make_figures.m.log');
fprintf('\nAnalysis started on %s\n\n', datetime);
clear;
close all;
addpath './matlab/classes';
addpath './matlab/functions';


%% Load data
data       = readtable('./data/intermediate/convolution-data.csv');
mle_output = load('./matlab/mat/mle.mat');
mle_output = mle_output.mle_output;
mle_output = mle_output(strcmp(mle_output.algorithm, 'gradient knitro'), :);


%% Create variables
% Prior
g          = @t_distribution;
dg         = @d_t_distribution;

g_normal = @(xx, bb) normpdf(xx, bb(1), bb(2));

% Select data for success rate
data = data(strcmp(data.Metric, 'session success rate'), :);

% Select MLE output with preferred specification of no weights.
row = (mle_output.use_weights == 0);
mle_output = mle_output(row, :);

% Choose beta for success rate with benchmark estimate
row = strcmp(mle_output.metric, 'session success rate');
beta = mle_output(row, 'beta');
beta = table2array(beta);

% Create useful variables
delta_hat = data.Delta;
sigma_n = data.StdErrorDelta;
sigma = data.sigma;

mean_sigma = mean(sigma);

% Sort variables by the signal
[delta_hat, sorted_index] = sort(delta_hat);
sigma = sigma(sorted_index);
sigma_n = sigma_n(sorted_index);

% Set normal distribution betas
beta_normal_shape = beta(1:2);
beta_normal_moments = [0, 0];
beta_normal_moments(1) = mean(delta_hat);
beta_normal_moments(2) = var(delta_hat) - mean(sigma_n .^ 2);
beta_normal_moments(2) = sqrt(beta_normal_moments(2));


%% model-fit-histogram.eps and model-fit-qq.eps
n_points = 400; % The density looks nicer if its precise.
grid_delta_hat = linspace(-0.3, 0.3, n_points);
% Calculate empirical density
theoretical_density = zeros(n_points, 1);
theoretical_cdf = zeros(n_points, 1);
for ii = 1:n_points
    theoretical_density(ii) = ...
        mean(Twee.L(grid_delta_hat(ii), sigma_n, beta , g));
    theoretical_cdf(ii) = ...
        mean(Twee.L_cdf(grid_delta_hat(ii), sigma_n, beta , g));    
end

% Plot histogram
figure();
set(gca,'FontSize',18);
hold on;

histogram(delta_hat, 'Normalization', 'pdf');
plot(grid_delta_hat, theoretical_density);

% Labels
xlabel({'Signal $\widehat{\delta}$'},'Interpreter','latex');
ylabel('Density');

% Save
print(gcf,'-depsc2','./output/figures/model-fit-histogram.eps');

% Plot qq
y  = zeros(length(delta_hat), 1);
for ii = 1:length(delta_hat)
    q = mean(Twee.L_cdf(delta_hat(ii), sigma_n, beta , g));
    y(ii) = quantile(delta_hat, q);
end

figure();
set(gca,'FontSize',18);
hold on;

scatter(delta_hat, y);
plot(grid_delta_hat, grid_delta_hat);

% Labels
xlabel('Theoretical quantile (fitted model)');
ylabel('Data quantile');

% Save
print(gcf,'-depsc2','./output/figures/model-fit-qq.eps');

%% Hill-plot.eps: Supplementary Materials

X = log(delta_hat(delta_hat>0));
X = sort(X,'descend');  %Sort values from largest to smallest

%Construct the Hill estimator as defined in the paper of Hausler 
%and Teugels AoS (1985)
% "On Asymptotic Normality of Hill's estimator for the exponent of 
% regular variation"

aux1 ...
  =cumsum(X)./((1:1:size(X,1))'); 
                    

Hillstat ...
  = aux1(1:end-1,1)-X(2:end,1);

figure()
hold on
plot(1:1:size(X,1)-1,1./Hillstat,'black','LineWidth',2);
plot([0 800],[beta(1,3) beta(1,3)],'--black','LineWidth',2)
hold off
xlabel('Number of top values used for estimation')
ylabel('Tail-index estimator')
ylim([0,5])
Hillaux=floor(.25*size(delta_hat,1));
xlim([0,Hillaux])

clear X aux1 Hillaux

print(gcf,'-depsc2','./output/figures/Hillplot.eps');

%% posterior-mean.eps: Posterior mean t-distribution
% Define the grid
grid_delta_hat = linspace(-0.2, 0.2, 100);
n_users = 20e6; 
typical_sigma_n = mean_sigma ./ sqrt(n_users);

y = Twee.mean_posterior(grid_delta_hat, typical_sigma_n, beta, g);


% Plot figure 1 for the paper
figure();
set(gca,'FontSize',18);

hold on;
plot(grid_delta_hat, y, 'black', 'LineWidth', 2.5);
plot(grid_delta_hat, 0.* grid_delta_hat, 'black');
plot(0 .* grid_delta_hat, grid_delta_hat, 'black');
hold off

% Labels
xlabel({'Signal $\widehat{\delta} _i$'},'Interpreter','latex');
ylabel({'Posterior mean $P_i(\widehat{\delta} _i, n_i)$'},'Interpreter','latex');

% Save
print(gcf,'-depsc2','./output/figures/posterior-mean.eps');

%% Supplementary Mateials: Posterior mean t-distribution vs. normal

y_normal_supp =  Twee.mean_posterior(grid_delta_hat, ...
                 typical_sigma_n, beta_normal_shape, g_normal);
             
figure();
set(gca,'FontSize',18);

hold on;
plot(grid_delta_hat, y, 'black', 'LineWidth', 2.5);
plot(grid_delta_hat, y_normal_supp, ':black', 'LineWidth', 2.5);
plot(grid_delta_hat, 0.* grid_delta_hat, 'black');
plot(0 .* grid_delta_hat, grid_delta_hat, 'black');
hold off

% Labels
xlabel({'Signal $\widehat{\delta} _i$'},'Interpreter','latex');
ylabel({'Posterior mean $P_i(\widehat{\delta} _i, n_i)$'},'Interpreter','latex');             

clear y_normal_supp

% Save
print(gcf,'-depsc2',...
      './output/figures/posterior-mean-supp-materials-t-vs-normal.eps');

%% Supplementary Materials: Posterior mean t-distribution (fat,edge,thin)

n_users_comparison = [20e6,10e6,5e6];

typical_sigma_n_comparison ...
                   = mean_sigma ./ sqrt(n_users_comparison);

% Calculate mean posterior for t-distribution vs. normal

y_t_dist_fat   = zeros(size(n_users_comparison,2),size(grid_delta_hat,2));

y_t_dist_edge  = zeros(size(n_users_comparison,2),size(grid_delta_hat,2));

y_t_dist_thin  = zeros(size(n_users_comparison,2),size(grid_delta_hat,2));

beta_thin      = [beta(1,1),beta(1,2),20];

beta_edge      = [beta(1,1),beta(1,2),3];

for i_y = 1: size(y_t_dist_fat,1)

y_t_dist_fat(i_y,:) ...
    = Twee.mean_posterior(grid_delta_hat, ...
                          typical_sigma_n_comparison(1,i_y), ...
                          beta, g);

y_t_dist_edge(i_y,:) ...
    = Twee.mean_posterior(grid_delta_hat, ...
                          typical_sigma_n_comparison(1,i_y), ...
                          beta_edge, g);
                                            
y_t_dist_thin(i_y,:) ...
    = Twee.mean_posterior(grid_delta_hat, ...
                          typical_sigma_n_comparison(1,i_y), ...
                          beta_thin, g);

end

%Plot figure for supplementary materials
figure();
set(gca,'FontSize',18);

hold on;
plot(grid_delta_hat, y_t_dist_fat(1,:), 'black', 'LineWidth', 2.5);
plot(grid_delta_hat, y_t_dist_fat(2,:), '--black', 'LineWidth', 2.5);
plot(grid_delta_hat, y_t_dist_fat(3,:), ':black', 'LineWidth', 2.5);
plot(grid_delta_hat, 0.* grid_delta_hat, 'black');
plot(0 .* grid_delta_hat, grid_delta_hat, 'black');
hold off

print(gcf,'-depsc2',...
     './output/figures/posterior-mean-supp-materials-fat.eps');

figure();
set(gca,'FontSize',18);

hold on;
plot(grid_delta_hat, y_t_dist_edge(1,:), 'black', 'LineWidth', 2.5);
plot(grid_delta_hat, y_t_dist_edge(2,:), '--black', 'LineWidth', 2.5);
plot(grid_delta_hat, y_t_dist_edge(3,:), ':black', 'LineWidth', 2.5);
plot(grid_delta_hat, 0.* grid_delta_hat, 'black');
plot(0 .* grid_delta_hat, grid_delta_hat, 'black');
hold off

print(gcf,'-depsc2',...
     './output/figures/posterior-mean-supp-materials-edge.eps');

figure();
set(gca,'FontSize',18);

hold on;
plot(grid_delta_hat, y_t_dist_thin(1,:), 'black', 'LineWidth', 2.5);
plot(grid_delta_hat, y_t_dist_thin(2,:), '--black', 'LineWidth', 2.5);
plot(grid_delta_hat, y_t_dist_thin(3,:), ':black', 'LineWidth', 2.5);
plot(grid_delta_hat, 0.* grid_delta_hat, 'black');
plot(0 .* grid_delta_hat, grid_delta_hat, 'black');
hold off

print(gcf,'-depsc2',...
     './output/figures/posterior-mean-supp-materials-thin.eps');

%% Constants related to the posterior mean
% sigma-n-typical
% representative standard error used in counterfactuals
% assumes the average sigma and 20 million users.
fid = fopen("./output/constants/sigma-n-typical.txt", "w");
    fprintf(fid, "%.3f", typical_sigma_n);
fclose(fid);

fid = fopen("./output/constants/marginally-significant-delta.txt", "w");
    fprintf(fid, "%.3f", 2*typical_sigma_n);
fclose(fid);

fid = fopen("./output/constants/outlier-delta.txt", "w");
    fprintf(fid, "%.3f", 4*typical_sigma_n);
fclose(fid);

fid = fopen("./output/constants/marginally-significant-p.txt", "w");
    p = Twee.mean_posterior(2*typical_sigma_n, typical_sigma_n, beta, g);
    fprintf(fid, "%.3f", p);
fclose(fid);

fid = fopen("./output/constants/outlier-p.txt", "w");
    p = Twee.mean_posterior(4*typical_sigma_n, typical_sigma_n, beta, g);
    fprintf(fid, "%.3f", p);
fclose(fid);

fid = fopen("./output/constants/typical-delta-star.txt", "w");
    typical_delta_star = Twee.mean_posterior_root(typical_sigma_n, beta, g);
    fprintf(fid, "%.3f", typical_delta_star);
fclose(fid);

fid = fopen("./output/constants/typical-t-star.txt", "w");
    typical_t_star = typical_delta_star / typical_sigma_n;
    fprintf(fid, "%.3f", typical_t_star);
fclose(fid);

fid = fopen("./output/constants/cost-to-rationalize-p-value.txt", "w");
    cost_to_rationalize = Twee.mean_posterior(1.96 * typical_sigma_n, typical_sigma_n, beta, g);
    fprintf(fid, "%.4f", cost_to_rationalize);
fclose(fid);


%% Pareto principle
top_percent = 0.02; % top 2% ideas.

% Theoretical sorted on delta
% This is commented because we will not use it in the paper. The result
% ends up being that over 90% of the gain comes from the top 2% and even 1%
% ideas. That actually makes sense because the prior is super concentrated
% around 0 and with fat tails. So it is really only outliers that matter.
% But let's do the more conservative and easier to understand data-driven
% approach for the paper.
% Define distribution of delta
% f = @(x) g(x, beta);
% F = @(x) Twee.robust_integral(f, [-inf, x, min(x, 0)]);
% 
% top_delta = fzero(@(x) F(x) - (1-top_percent), 0);
% 
% % Define probability idea gets implemented under 5% significance rule
% probability_implemented = @(x) ...
%     normcdf(x / typical_sigma_n - 1.96);
% 
% % Calculate total payoff of 5% rule
% waypoints = ...
%     [0, inf, ...
%     typical_sigma_n, 2*typical_sigma_n, 3*typical_sigma_n, 5*typical_sigma_n, 10*typical_sigma_n, ...
%     beta(2), 2*beta(2), 3*beta(2), 5*beta(2), 10*beta(2)];
% waypoints = [waypoints, -waypoints];
% total_payoff = Twee.robust_integral(...
%     @(x) x .* probability_implemented(x) .* f(x), waypoints);
% waypoints = max(waypoints, top_delta);
% top_payoff = Twee.robust_integral(...
%     @(x) x .* probability_implemented(x) .* f(x), waypoints);
% 
% pareto_ratio_theory = 100 * top_payoff / total_payoff;


% Data driven
p_hat = Twee.mean_posterior(delta_hat, sigma_n, beta, g);

policy_list = {"p-value", "bayesian"};
eval_type_list = {"naive", "bayesian"};
gains_table = table();
threshold_top = quantile(delta_hat, 1 - top_percent);
is_top = delta_hat > threshold_top;
current_row = 0;


for ii = 1:numel(policy_list)
    policy = policy_list{ii};
    if strcmp(policy, "p-value")
        is_implemented = delta_hat > 1.96 * sigma_n;
    elseif strcmp(policy, "bayesian")
        is_implemented = p_hat> 0;
    end
    for jj = 1:numel(eval_type_list)
        eval_type = eval_type_list{jj};
        if strcmp(eval_type, "naive")
            gain = delta_hat;
        elseif strcmp(eval_type, "bayesian")
            gain = p_hat;
        end
        
        current_row = current_row + 1;
        gains_table(current_row, {'spec'}) = {strcat(policy, " ", eval_type)};
        gains_table(current_row, {'total_gain'}) = {sum(gain .* is_implemented)};
        gains_table(current_row, {'top_gain'}) = {sum(gain .* is_implemented .* is_top)};
    end
end
gains_table.pareto_ratio = 100 * gains_table.top_gain ./ gains_table.total_gain;

% Constants
gains_optimal = ...
    gains_table.total_gain(strcmp("bayesian bayesian", gains_table.spec));
gains_p_value = ...
    gains_table.total_gain(strcmp("p-value bayesian", gains_table.spec));
optimal_vs_p_value_gain = 100 * (gains_optimal - gains_p_value) / gains_p_value;
pareto_ratio_p_value = ...
    gains_table.pareto_ratio(strcmp("p-value bayesian", gains_table.spec));

fid = fopen("./output/constants/pareto-ratio-p-value.txt", "w");
    fprintf(fid, "%.1f", pareto_ratio_p_value);
fclose(fid);

fid = fopen("./output/constants/historical-gain-optimal.txt", "w");
    fprintf(fid, "%.1f", gains_optimal);
fclose(fid);

fid = fopen("./output/constants/optimal-vs-p-value-gain.txt", "w");
    fprintf(fid, "%.2f", optimal_vs_p_value_gain);
fclose(fid);

%% production-function-supplementary material: 

% Graph parameters
n_left = 0;
n_right = 40e6;
n_points_alpha = 3;
values_alpha = [beta(1,3),3,20];
n_points_production_function = 100;

% Set variables
grid_n = linspace(n_left, n_right, n_points_production_function);
y_t_supp = zeros(n_points_production_function, n_points_alpha);

for jj = 1:n_points_alpha
    f_infinity = Twee.f(inf, mean_sigma, beta + (jj - 1) * [0 0 step_alpha], g);
    for ii = 1:n_points_production_function
    n = grid_n(ii);
    y_t_supp(ii, jj) = ...
        Twee.f(n, mean_sigma, ...
              [beta(1,1),beta(1,2),values_alpha(1,jj)], g) / f_infinity;
    end
end

figure();
    plot(grid_n / 1e6, y_t_supp(:,1),'black');
    set(gca, 'FontSize', 18);
    grid off;
    xlabel('Size of the experiment $n$ (millions)', 'interpreter', 'latex');
    ylabel('Production function $f(n)/f(\infty)$', 'interpreter', 'latex');
    set(gca, 'YTickMode','manual');
    set(gca, 'YTickLabel', ...
        num2str(100.*get(gca, 'YTick')', ...
        '%g%%'));
    
print(gcf,'-depsc2',...
     './output/figures/production-function-supp-materials-fat.eps');    
    
figure();
    plot(grid_n / 1e6, y_t_supp(:,2),'black');
    set(gca, 'FontSize', 18);
    grid off;
    xlabel('Size of the experiment $n$ (millions)', 'interpreter', 'latex');
    ylabel('Production function $f(n)/f(\infty)$', 'interpreter', 'latex');  
    set(gca, 'YTickMode','manual');
    set(gca, 'YTickLabel', ...
        num2str(100.*get(gca, 'YTick')', ...
        '%g%%'));

print(gcf,'-depsc2',...
     './output/figures/production-function-supp-materials-edge.eps');       
        
 figure();
    plot(grid_n / 1e6, y_t_supp(:,3),'black');
    set(gca, 'FontSize', 18);
    grid off;
    xlabel('Size of the experiment $n$ (millions)', 'interpreter', 'latex');
    ylabel('Production function $f(n)/f(\infty)$', 'interpreter', 'latex'); 
    set(gca, 'YTickMode','manual');
    set(gca, 'YTickLabel', ...
        num2str(100.*get(gca, 'YTick')', ...
        '%g%%'));

print(gcf,'-depsc2',...
     './output/figures/production-function-supp-materials-thin.eps'); 
 
%% production-function-t.eps and production-function-normal.eps
% Graph parameters
n_left = 0;
n_right = 40e6;
n_points_alpha = 4;
step_alpha = 1;
n_points_production_function = 100;

% Set variables
grid_n = linspace(n_left, n_right, n_points_production_function);
y_t = zeros(n_points_production_function, n_points_alpha);
y_normal = zeros(n_points_production_function, 1);

% Calculate production function
% t distribution
for jj = 1:n_points_alpha
    f_infinity = Twee.f(inf, mean_sigma, beta + (jj - 1) * [0 0 step_alpha], g);
    for ii = 1:n_points_production_function
    n = grid_n(ii);
    y_t(ii, jj) = ...
        Twee.f(n, mean_sigma, beta + (jj - 1) * [0 0 step_alpha], g) / f_infinity;
    end
end

% normal distribution
f_infinity = Twee.f(inf, mean_sigma, beta_normal_shape, g_normal);
for ii = 1:n_points_production_function
    n = grid_n(ii);
    y_normal(ii) = ...
        Twee.f(n, mean_sigma, beta_normal_shape, g_normal) / f_infinity;
end

% Plot t distribution
figure();
    plot(grid_n / 1e6, y_t);
    set(gca, 'FontSize', 18);
    grid off;
    xlabel('Size of the experiment $n$ (millions)', 'interpreter', 'latex');
    ylabel('Production function $f(n)/f(\infty)$', 'interpreter', 'latex');
    
    % Label curves
    for jj = 1:n_points_alpha
        ii = floor(3 / 4 * n_points_production_function);
        n = grid_n(ii) / 1e6;
        y = y_t(ii, jj) + ...
            -0.030 * max(y_t(:)); % vertical offset
        ss = num2str(beta(3) + (jj - 1) * step_alpha, 3);
        ss = ['$\alpha = ', ss, '$'];
        text(n, y, ss, 'interpreter', 'latex', 'FontSize', 18);
    end
    % Axis in percent
    set(gca, 'YTickMode','manual');
    set(gca, 'YTickLabel', ...
        num2str(100.*get(gca, 'YTick')', ...
        '%g%%'));
    
% Save
print(gcf,'-depsc2','./output/figures/production-function-t.eps');

% Plot normal distribution
figure();
    plot(grid_n / 1e6, y_normal);
    set(gca,'FontSize',18);
    grid off;
    % title('Normal distribution', 'interpreter', 'latex');
    xlabel('Size of the Experiment $n$ (millions)', 'interpreter', 'latex');
    ylabel('Production function $f(n)/f(\infty)$', 'interpreter', 'latex');
    % Axis in percent
    set(gca, 'YTickMode','manual')
    set(gca, 'YTickLabel', ...
        num2str(100.*get(gca, 'YTick')', ...
        '%g%%'));
    
% Save
print(gcf,'-depsc2','./output/figures/production-function-normal.eps');    


%% Constants: counterfactuals of increasing number of ideas or data
n_users = 20e6;
I = length(delta_hat);
n_users_total = I * n_users;

current_f = Twee.f(n_users, mean_sigma, beta, g);

% ten percent more ideas
leaner_10_f = 1.1 * Twee.f(n_users / 1.1, mean_sigma, beta, g);
leaner_10_gain = 100 * (leaner_10_f - current_f) / current_f;
fid = fopen("./output/constants/lean-gain-10.txt", "w");
    fprintf(fid, "%.2f", leaner_10_gain);
fclose(fid);

% twenty percent more ideas
leaner_20_f = 1.2 * Twee.f(n_users / 1.2, mean_sigma, beta, g);
leaner_20_gain = 100 * (leaner_20_f - current_f) / current_f;
fid = fopen("./output/constants/lean-gain-20.txt", "w");
    fprintf(fid, "%.2f", leaner_20_gain);
fclose(fid);

% ten percent more data
current_F = I * current_f;
more_data_F = I * Twee.f(n_users * 1.1, mean_sigma, beta, g);
more_data_gain_F = more_data_F - current_F;
more_data_gain_pc = 100 * more_data_gain_F / current_F;
fid = fopen("./output/constants/more-data-gain-F.txt", "w");
    fprintf(fid, "%.2f", more_data_gain_F);
fclose(fid);
fid = fopen("./output/constants/more-data-gain-pc.txt", "w");
    fprintf(fid, "%.2f", more_data_gain_pc);
fclose(fid);


%% Constants (others)
% alpha
% benchmark estimate of alpha
fid = fopen("./output/constants/alpha.txt", "w");
    fprintf(fid, "%.2f", beta(3));
fclose(fid);

% alternative alpha
alternative_alpha = beta(3) + 1;
fid = fopen("./output/constants/alternative-alpha.txt", "w");
    fprintf(fid, "%.2f", alternative_alpha);
fclose(fid);

% Decrease with alternative alpha
loss_alternative_alpha = 1 - ...
    Twee.f(inf, mean_sigma, [beta(1:2) alternative_alpha], g) / ...
    Twee.f(inf, mean_sigma, beta, g);
loss_alternative_alpha = 100 * loss_alternative_alpha;
fid = fopen("./output/constants/loss-alternative-alpha.txt", "w");
    fprintf(fid, "%.2f", loss_alternative_alpha);
fclose(fid);

% Probability of an improvement better than 10%
prob_large_gain = integral(@(d) g(d, beta), 10, Inf);
fid = fopen("./output/constants/probability-more-than-10pc-gain.txt", "w");
    fprintf(fid, "%.2e", prob_large_gain);
fclose(fid);

% Value f of a typical experiment
typical_f = Twee.f(20e6, mean_sigma, beta, g);
fid = fopen("./output/constants/typical-f.txt", "w");
    fprintf(fid, "%.2e", typical_f);
fclose(fid);


%% mle-results-mean-scale.eps and mle-results-tail.eps
%% MUST COME LAST BECAUSE THIS CODE MESSES UP ALL VARIABLES!
% Convert metric names to names in the paper
mle_output.index = zeros(height(mle_output), 1);
for ii = 1:height(mle_output)
    switch mle_output.metric{ii}
        case 'session success rate'
            mle_output.metric{ii} = 'SR0';
            mle_output.index(ii) = 1;
        case 'page click rate'
            mle_output.metric{ii} = 'SR1';
            mle_output.index(ii) = 2;            
        case 'quickback rate'
            mle_output.metric{ii} = 'SR2';
            mle_output.index(ii) = 3;            
        case 'time to success'
            mle_output.metric{ii} = 'SR3';
            mle_output.index(ii) = 4;            
        case 'queries per user'
            mle_output.metric{ii} = 'LR1';
            mle_output.index(ii) = 5;            
        case 'sessions per user'
            mle_output.metric{ii} = 'LR2';
            mle_output.index(ii) = 6;
    end
end

% Drop unused metrics.
mle_output(mle_output.index == 0, :) = [];

% Order by metrics
mle_output = sortrows(mle_output, 'index');

% Plot tail coefficients
figure();
set(gca, 'FontSize', 18);
hold on;

for ii = 1:height(mle_output)
    beta = mle_output(ii, 'beta');
    beta = table2array(beta);
    omega = mle_output(ii, 'variance_matrix');
    omega = omega{1, 1};
    omega = omega{1};

    alpha = beta(3);
    c = 1.96 * sqrt(omega(3, 3));
    ub = alpha + c;
    lb = alpha - c;
    
    % Plot confidence intervals
    plot([ii ii], [lb ub], '-b');
    % Plot point estimates
    plot(ii, alpha, 'o','MarkerFaceColor','b', 'MarkerEdgeColor', 'b');
    % Text for estimate
    text(ii, alpha, ['  ' num2str(round(alpha,2))], 'FontSize', 14);
end

xlim([0, height(mle_output) + 1]);
ylim([0, 4.5]);
set(gca, 'XTick', 1:6, 'XTickLabel', mle_output.metric)

% Plot horizontal line at alpha = 3
current_axis = gca;
plot(current_axis.XLim, [3 3], ':b');

xlabel('Metric')
ylabel('Tail coefficient \alpha')

% Save
print(gcf,'-depsc2','./output/figures/mle-results-tail.eps')

% Plot mean and scale
figure();
set(gca,'FontSize',18);
ylim([-2e-3, 15e-3]);
hold on;

for ii = 1:height(mle_output)
    beta = mle_output(ii, 'beta');
    beta = table2array(beta);
    omega = mle_output(ii, 'variance_matrix');
    omega = omega{1, 1};
    omega = omega{1};

    M = beta(1);
    c = 1.96 * sqrt(omega(1, 1));
    M_ub = M + c;
    M_lb = M - c;
    s = beta(2);
    c = 1.96 * sqrt(omega(2, 2));
    s_ub = s + c;
    s_lb = s - c;
    
    % Plot confidence intervals
    plot([M M], [s_lb s_ub], 'b');
    plot([M_lb M_ub], [s s], 'b');
    % Plot point estimates
    plot(M, s, 'o','MarkerFaceColor','b', 'MarkerEdgeColor', 'b');
    % Text for estimate
    if ~strcmp(mle_output.metric{ii}, 'LR2')
        text(M + 0.5e-3, s + 0.5e-3, mle_output.metric{ii}, 'FontSize', 14);
    else
        text(M - 0.5e-3, s + 0.5e-3, mle_output.metric{ii}, 'FontSize', 14);
    end
end

% Plot lines for zero
current_axis = gca;
plot(current_axis.XLim, [0 0], ':b');
plot([0 0], current_axis.YLim, ':b');

xlabel('Mean $M$', 'Interpreter', 'LaTex')
ylabel('Scale $s$', 'Interpreter', 'LaTex')

% Save
print(gcf,'-depsc2','./output/figures/mle-results-mean-scale.eps')


%% End
fprintf('\nAnalysis ended on %s\n\n', datetime);
diary off;