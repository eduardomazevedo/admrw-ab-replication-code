function [f] = f_gaussian(n, sigma, beta)
%F_GAUSSIAN Production function of the gaussian gaussian model
%   f_gaussian(n, sigma, beta) gives the production function with n
%   samples, experimental noise standard error sigma, and gaussian parameter vector
%   beta with mean and standard error.
m = beta(1);
s = beta(2);

theta = @(nn) s .^ 2 ./ (s.^2 + sigma.^2 ./ nn);
conditional_var = @(nn) theta(nn) .* sigma .^ 2 ./ nn;

mp = @(nn) 1 ./ 2 ./ nn .* ...
    normpdf(m ./ theta(nn), 0, sqrt(s.^2 + sigma.^2 ./ nn)) .* ...
    conditional_var(nn);

f = integral(mp, 0, n);

end

