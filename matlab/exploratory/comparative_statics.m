%% Exploratory analysis of comparative statics

close all;

% Graph
beta = beta_normal_shape .* [1 1.7];
n_points = 100;

n_grid = linspace(0, 4e7, n_points);

for ii = 1:n_points
    f_infinity = Twee.f(Inf, mean_sigma, beta, g_normal);
    n = n_grid(ii);
    y(ii) = f_gaussian(n, mean_sigma, beta) / f_infinity;
end

plot(n_grid, y);


%% Average product maximization
mm = 1;
beta = beta_normal_shape .* [1 6];
ap = @(nn) f_gaussian(nn, mean_sigma, [beta(1), beta(2) * mm]) ./ nn - 1e-3 ./ nn;

efficient_scale = fminsearch(@(nn) -ap(nn), 20e6)


%% Average product maximization
n_points = 100;

m_grid = linspace(0, 20, n_points);
s_grid = linspace(0.01, 4, n_points);

z = zeros(n_points, n_points);
for ii = 1:n_points
    for jj = 1:n_points
        beta = [m_grid(ii), s_grid(jj)];
        ap = @(nn) f_gaussian(nn, mean_sigma, [beta(1), beta(2)]) ./ nn;
        efficient_scale = fminsearch(@(nn) -ap(nn), 2e4);
        z(ii, jj) = efficient_scale;
    end
end

for ii = 1:n_points
    decreasing(ii) = issorted(-z(ii, :));
end