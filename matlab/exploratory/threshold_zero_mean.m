%% Start
clear;
addpath 'classes'
addpath 'functions'
load('./mat/mle.mat');
close all;

%% Constants
n_points = 500;

%% Define basics
beta = mle_output(7, :).beta;
clear mle_output;

g = @t_distribution;
f = @(nn, bb) Twee.f(nn, 100, bb, g);
z_star = @(nn, bb) Twee.mean_posterior_root(100 / sqrt(nn), bb, g);

% Zero mean
beta(1) = 0;

%% Plots
n_grid = exp(linspace(log(1e3), log(1e6), n_points))';
v_grid = linspace(1.5, 5, 5);
z_grid = zeros(n_points, 5);
f_grid = zeros(n_points, 5);

for i = 1:n_points
    for j = 1:5
        b = beta;
        b(3) = v_grid(j);
        z_grid(i, j) = z_star(n_grid(i), b);
        f_grid(i, j) =      f(n_grid(i), b);
    end
end

% Plot times n
figure();
plot(n_grid, z_grid .* n_grid);
legend(num2str(v_grid));
title('z*n');

% Plot times square root of n
figure();
plot(n_grid, z_grid .* sqrt(n_grid));
legend(num2str(v_grid));
title('z*sqrt(n)');

% Log log plot of z
figure();
loglog(n_grid, z_grid);
legend(num2str(v_grid));
title('z loglog');

% Log log plot of f
figure();
loglog(n_grid, f_grid);
legend(num2str(v_grid));
title('f loglog');

% Regular plot of f
figure();
plot(n_grid, f_grid);
legend(num2str(v_grid));
title('f');

% Regular plot of f
figure();
plot(n_grid, f_grid(:, 2:5));
legend(num2str(v_grid(2:5)));
title('f high betas');

% Marginal product
df_grid = zeros(n_points, 5);
for j = 1:5
    df_grid(:, j) = gradient(f_grid(:, j), n_grid);
end

% Log log plot of df
figure();
loglog(n_grid, df_grid);
legend(num2str(v_grid));
title('df loglog');

% Regular plot of df
figure();
plot(n_grid, df_grid);
legend(num2str(v_grid));
title('df');

% Regular plot of f
figure();
plot(n_grid, df_grid(:, 2:5));
legend(num2str(v_grid(2:5)));
title('df high betas');
