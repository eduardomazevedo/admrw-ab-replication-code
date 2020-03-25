%% Start
addpath classes;
addpath functions;
clear;
close all;

% Print date
fprintf('\nAnalysis started on %s\n\n', datetime);

% Set parameters to SSR
beta = [-0.00044808, 0.0044135, 1.7189];
sigma0 = 100;
% beta = beta + 4;

% Set parameters to sessions per user
% beta = [-0.00020467 8.7364e-05 4.0425];
% sigma0 = 300;

g  = @t_distribution;
dg = @d_t_distribution;


%% Likelihood
for ii = 1:3
    % Grid for plotting likelihood
    if ii == 1
        sigma = sigma0;
    elseif ii == 2
        sigma = beta(2);
    elseif ii == 3
        sigma = beta(2) / 10;
    end
    z_grid = linspace(-1, 1) * 6 * sigma;

    % Plot Likelihood
    figure();
    plot(z_grid, Twee.L(z_grid, sigma, beta, g));

    % Plot log likelihood
    figure();
    plot(z_grid, Twee.l(z_grid, sigma, beta, g));

    % Plot cumulative likelihood
    figure();
    plot(z_grid, Twee.L_cdf(z_grid, sigma, beta, g));
end;


%% Posterior mean
for ii = 1:3
    % Grid for plotting
    if ii == 1
        sigma = sigma0;
    elseif ii == 2
        sigma = beta(2);
    elseif ii == 3
        sigma = beta(2) / 10;
    end
    z_grid = linspace(-1, 1) * 6 * sigma;

    % Plot posterior
    figure();
    plot(z_grid, Twee.mean_posterior(z_grid, sigma, beta, g));
end;


%% Posterior mean root
for ii = 1:6
    % Grid for plotting
    min_n = 10^(ii-1);
    max_n = 10^(ii+1);
    n_grid = linspace(min_n, max_n);

    % Plot posterior
    figure();
    tic;
    plot(n_grid, Twee.mean_posterior_root(sigma0 ./ sqrt(n_grid), beta, g));
    toc
end;


%% Maximum of f versus posterior mean root
close all;
for power = 1:7
    n = 10^power;

    z_bar_root = Twee.mean_posterior_root(sigma0 / sqrt(n), beta, g);
    f_z_bar_root = Twee.f(n, sigma0, beta, g, z_bar_root);
    [best_f, z_bar_optimization] = Twee.f(n, sigma0, beta, g);

    z_grid = linspace(0, 2 * z_bar_root);
    z_grid = sort([z_grid z_bar_root z_bar_optimization]);

    % Graphs calculated point by point
    f_graph_one_by_one = 0 .* z_grid;
    posterior_graph_one_by_one = 0 .* z_grid;
    for ii = 1:length(z_grid)
        f_graph_one_by_one(ii) = Twee.f(n, sigma0, beta, g, z_grid(ii));
        posterior_graph_one_by_one(ii) = Twee.mean_posterior(z_grid(ii), sigma0 / sqrt(n), beta, g);
    end

    % Plot f
    figure();
    hold on;
    plot(z_grid, max(Twee.f(n, sigma0, beta, g, z_grid), 0));
    axis = gca();
    line([z_bar_root z_bar_root], [axis.YLim(1), axis.YLim(2)]);
    line([z_bar_optimization z_bar_optimization], [axis.YLim(1), axis.YLim(2)], 'Color', 'red');
    line(z_grid, z_grid.*0 + best_f, 'Color', 'green');
    line(z_grid, z_grid.*0 + f_z_bar_root, 'Color', 'red');
    line(z_grid, max(0, f_graph_one_by_one), 'Color', 'red');

    figure();
    hold on;
    plot(z_grid, [Twee.mean_posterior(z_grid, sigma0 / sqrt(n), beta, g); z_grid .*0]);
    axis = gca();
    line([z_bar_root z_bar_root], [axis.YLim(1), axis.YLim(2)]);
    line([z_bar_optimization z_bar_optimization], [axis.YLim(1), axis.YLim(2)], 'Color', 'red');
    line(z_grid, posterior_graph_one_by_one, 'Color', 'red');

    hold off;
end


%% Production function
for ii = 1:6
    % Grid for plotting
    min_n = 10^(ii-1);
    max_n = 10^(ii+1);
    n_grid = linspace(min_n, max_n);
    
    % Calculate f
    yy = n_grid .* 0;
    for jj = 1:length(n_grid)
        yy(jj) = Twee.f(n_grid(jj), sigma0, beta, g);
    end

    % Plot posterior
    figure();
    plot(n_grid, yy);
end;


%% Production function implementation details
Twee.f(1e6, sigma0, beta, g)

Twee.f(Inf, sigma0, beta, g)

try
    Twee.f([1e6 1e7], sigma0, beta, g);
catch
    display('correctly returned an error with non scalar input');
end

try
    Twee.f(1e6, [sigma0, sigma0+1], beta, g);
catch
    display('correctly returned an error with non scalar input');
end


%% End
fprintf('\nAnalysis started on %s\n\n', datetime);
close all;