classdef Twee
    %Tweedie collection of functions to do Tweedie's formula related tasks.
    %
    methods(Static = true)
        function I = robust_integral(integrand, waypoints)
            % ROBUST_INTEGRAL Robust integral
            % robust_integral(integrand, waypoints) evaluates the integral
            % adding integrals for each waypoint. This is slower but
            % numerically more robust than the built-in MATLAB function.
            % The range of integration is the range of waypoints.
            waypoints = sort(waypoints);
            output_dimensions = size(integrand(0));
            I = zeros(output_dimensions);
            for ii = 1:length(waypoints) - 1
                I = I + ...
                    integral(integrand, ...
                    waypoints(ii), waypoints(ii+1), 'ArrayValued', true);
            end
        end
        
        function L = L(z, sigma, beta, g)
            % L Likelihood.
            %   L(z, sigma, beta, g) calculates the likelihood of signal z,
            %   given signal standard deviation sigma, parameters beta and
            %   delta distribution g(delta, beta).
            integrand = @(delta) ...
                normpdf(z, delta, sigma) .* g(delta, beta);
            % Choose waypoints.
            zz = z(:)';
            waypoints = [ ...
                -Inf, Inf, ...
                beta(1), ...
                beta(1) - 10*beta(2), beta(1) + 10*beta(2), ...
                beta(1) - 10*min(sigma(:)), beta(1) + 10*min(sigma(:)), ...
                beta(1) - 10*max(sigma(:)), beta(1) + 10*max(sigma(:)), ...
                zz, ...
                zz - 10*beta(2), zz + 10*beta(2), ...
                zz - 10*min(sigma(:)), zz + 10*min(sigma(:)), ...
                zz - 10*max(sigma(:)), zz + 10*max(sigma(:)), ...
                0, ...
                - 10*beta(2), 10*beta(2), ...
                - 10*min(sigma(:)), 10*min(sigma(:)), ...
                - 10*max(sigma(:)), 10*max(sigma(:)), ...
                ];
            % Calculate integral
            L = Twee.robust_integral(integrand, waypoints);
        end

        function l = l(z, sigma, beta, g)
            % L log likelihood.
            %   L(z, sigma, beta, g) calculates the log-likelihood of signal z,
            %   given signal standard deviation sigma, parameters beta and
            %   delta distribution g(delta, beta).
            l = log(Twee.L(z, sigma, beta, g));
        end

        function L_cdf = L_cdf(z, sigma, beta, g)
            % L_CDF Likelihood CDF
            % L_cdf(z, sigma, beta, g) is the cumulative distribution of
            % the likelihood L(z, sigma, beta, g).
            integrand = @(delta) ...
                normcdf(z, delta, sigma) .* g(delta, beta);
            % Choose waypoints.
            zz = z(:)';
            waypoints = [ ...
                -Inf, Inf, ...
                beta(1), ...
                beta(1) - 10*beta(2), beta(1) + 10*beta(2), ...
                beta(1) - 10*min(sigma(:)), beta(1) + 10*min(sigma(:)), ...
                beta(1) - 10*max(sigma(:)), beta(1) + 10*max(sigma(:)), ...
                zz, ...
                zz - 10*beta(2), zz + 10*beta(2), ...
                zz - 10*min(sigma(:)), zz + 10*min(sigma(:)), ...
                zz - 10*max(sigma(:)), zz + 10*max(sigma(:)), ...
                0, ...
                - 10*beta(2), 10*beta(2), ...
                - 10*min(sigma(:)), 10*min(sigma(:)), ...
                - 10*max(sigma(:)), 10*max(sigma(:)), ...
                ];
            L_cdf = Twee.robust_integral(integrand, waypoints);
        end

        function delta_bar = mean_posterior(z, sigma, beta, g)
            % MEAN_POSTERIOR posterior mean of the true effect.
            %   MEAN_POSTERIOR(z, sigma, beta, g) gives the Bayesian estimate of the
            %   mean of delta given signal z, standard deviation of the
            %   signal sigma, parameters beta and distribution g(delta, beta) of delta.
            integrand = @(delta) ...
                delta .* ...
                normpdf(z, delta, sigma) .* g(delta, beta);
            % Choose waypoints.
            zz = z(:)';
            waypoints = [ ...
                -Inf, Inf, ...
                beta(1), ...
                beta(1) - 10*beta(2), beta(1) + 10*beta(2), ...
                beta(1) - 10*min(sigma(:)), beta(1) + 10*min(sigma(:)), ...
                beta(1) - 10*max(sigma(:)), beta(1) + 10*max(sigma(:)), ...
                zz, ...
                zz - 10*beta(2), zz + 10*beta(2), ...
                zz - 10*min(sigma(:)), zz + 10*min(sigma(:)), ...
                zz - 10*max(sigma(:)), zz + 10*max(sigma(:)), ...
                0, ...
                - 10*beta(2), 10*beta(2), ...
                - 10*min(sigma(:)), 10*min(sigma(:)), ...
                - 10*max(sigma(:)), 10*max(sigma(:)), ...
                ];
            delta_bar = ...
                Twee.robust_integral(integrand, waypoints) ./ ...
                Twee.L(z, sigma, beta, g);
        end

        function z = mean_posterior_root(sigma, beta, g)
            % MEAN_POSTERIOR_ROOT root of the mean posterior map.
            %   MEAN_POSTERIOR_ROOT(sigma, beta, g) returns the signal that implies in the mean posterior being zero.
            if isscalar(sigma)
                tolerance = min([beta(2), sigma, 1e-6]) * 1e-3;
                f = @(z) Twee.mean_posterior(z, sigma, beta, g);
%                 [z, error_f] = fzero(f, -beta(1)); Old version without
%                 the bisection method. Bisection is better for the cases
%                 where the posterior goes up too rapidly for newton's
%                 method.
                ub = sigma;
                while(f(ub) < 0)
                    ub = 10*ub;
                end
                lb = -sigma;
                while(f(lb) > 0)
                    lb = 10 * lb;
                end
                
                addpath('matlab/vendor/bisection');
                options = optimset('TolX', tolerance);
                [z, error_f, exit_flag] = bisection(f, lb, ub, 0, options);
                if exit_flag ~= 1
                    warning('Posterior inverse error!');
                end
            else
                f = @(x) Twee.mean_posterior_root(x, beta, g);
                z = arrayfun(f, sigma);
            end
        end

        function [beta_out, l, flag, output] = fit_g_no_gradient(z_data, sigma_data, beta_initial, g, max_iterations)
            % FIT_G_NO_GRADIENT estimates distribution of true effects.
            %   FIT_G_NO_GRADIENT(z_data, sigma_data, beta_initial, g) does maximum
            %   likelihood estimation of the parameter vector beta given data
            %   on z and sigma, initial value of beta, and functional form
            %   g(delta, beta).
            obj = @(beta) -mean(Twee.l(z_data, sigma_data, beta, g));
            options = optimset('Display', 'iter', 'MaxIter', max_iterations);
            [beta_out, l, flag, output] = ...
                fminsearch(obj, beta_initial, options);
            l = -l;
        end

        function [beta_out, ...
                  l       , ...
                  flag    , ...
                  output  , ...
                  MLasyvar] = fit_g(...
                    z_data, sigma_data, beta_initial, ...
                    g, dg, dgg, ...
                    max_function_evaluations, use_knitro, mle_weights)
            % FIT_G estimates distribution of true effects.
            %
            %   FIT_G(z_data, sigma_data, beta_initial, g, ...) does maximum
            %   likelihood estimation of the parameter vector beta given data
            %   on z and sigma, initial value of beta, and functional form
            %   g(delta, beta). 
            %
            %   Additional parameters: gradient dg and Hessian ddg,
            %   max_function_evaluations to be passed to optimizer,
            %   use_knitro flag to choose between knitro and MATLAB
            %   optimizer, mle_weights to use weighted observations.
            %
            %   Output: beta_out estimate, l value of log likelihood, flag
            %   returned by optimizer, output TODO, MLasyvar TODO
            
            % Define objective function.
            function [f, df, ddf, Omega] = obj(beta)
                                                                                                  
                LL      =  Twee.L(z_data', sigma_data', beta, g); 
                           %likelihood for each value of z
                           
                f       = -mean(mle_weights' .* log(LL), 2);
                           % Negative log-likelihood objective function
                           % (mean over all coordinates)
                
                df1     =  Twee.L(z_data', sigma_data', beta, dg);
                           % first-derivative of the likelihood with
                           % respect to beta
                
                df      = -mean(mle_weights' .* df1./LL, 2);
                           %sample score
                           
                if nargout > 2
                    
                d       =  size(beta, 2); 
                           %number of parameters to estimate           
                                           
                auxdf1  =  reshape(mle_weights'.*df1./LL, ...
                                   [d, 1, size(z_data, 1)]);
                           %reshape the individual scores
                
                Omega   = mean(bsxfun(@times, auxdf1, permute(auxdf1, [2, 1, 3])), 3);          
                           %sample variance of the score 
                
                dggaux  = @(x, beta) reshape(dgg(x, beta), [d^2, 1]);
                           %vectorization of the hessian of g
                
                df2     =  Twee.L(z_data', sigma_data', beta, dggaux);
                           %vectorized second-derivative of L
                
                ddf     = -reshape(mean(mle_weights' .* df2./LL, 2), ...
                           [size(beta, 2),size(beta, 2)]) ...
                          + Omega;
                           %sample information matrix 
                end
            end

            if ~exist('knitromatlab', 'file') || ~use_knitro
                options = optimoptions( ...
                    'fmincon','SpecifyObjectiveGradient',true, 'MaxFunctionEvaluations', max_function_evaluations,...
                    'Display','iter-detailed','Algorithm','trust-region-reflective', 'CheckGradients', false);
                [beta_out, l, flag, output] = fmincon(@obj, beta_initial, [], [], [], [], [-Inf, 0, 0], [], [], options);
                l = -l;
            else
                options = optimset('GradObj', 'on', 'Display', ...
                                   'iter-detailed');
                [beta_out, l, flag, output] = ...
                    knitromatlab(@obj, beta_initial, [], [], [], [], [-Inf, 0, 0], [], ...
                                 [], [], options);
                l = -l;
            end
            
            if nargout > 4
                
                [~, ~, ddf, Omega] = obj(beta_out);
                
                MLasyvar    = (ddf^(-1)) * (Omega) * (ddf^(-1)) / size(z_data, 1);
                
            end
        end

        function [f, z_bar] = f(n, sigma_individual, beta, g, varargin)
            % F production function
            %   F(n, sigma_individual, beta, g, [optional z_bar]) outputs the production
            %   function evaluated at n observations, individual experiment
            %   standard deviation sigma_individual, parameters beta and
            %   g(delta, beta) the true distribution of delta. The default
            %   is for the optimal production function. But if the user
            %   inputs an extra argument for a threshold z_bar it outputs
            %   that production function.
            
            % Check if arguments are scalar
            if ~isscalar(n) && isscalar(sigma_individual)
                error('f takes scalar n and sigma_individual arguments.');
            end
            
            % If n = 0 return 0, and if n = infinity return limit.
            if n == 0 || sigma_individual == Inf
                f = 0;
                return;
            elseif isinf(n)
                if isempty(varargin)
                    z_bar = 0;
                else
                    z_bar = varargin{1};
                end
                
                integrand = @(delta) delta .* g(delta, beta);
                waypoints = [0, Inf,...
                    beta(2), 10 * beta(2) ...
                    ];
                f = Twee.robust_integral(integrand, waypoints) - max(beta(1), 0);
                return;
            end
            
            % Define variables.
            m = beta(1);
            sigma = sigma_individual ./ sqrt(n);
            
            % If no extra threshold arguments, then calculate the efficient
            % production function.
            if isempty(varargin)
%                 z_bar = Twee.mean_posterior_root(sigma, beta, g);
                f_aux = @(z) -Twee.f(n, sigma_individual, beta, g, z);
                z_guess = Twee.mean_posterior_root(sigma, beta, g);
                [z_bar, fval] = fminsearch(f_aux, z_guess);
            elseif length(varargin) == 1
                z_bar = varargin{1};
            else
                error('Too many arguments!');
            end
            
            % Calculate f at z_bar
            % Define integral
            integrand = @(delta) ...
                delta .* ...
                (1 - normcdf(z_bar, delta, sigma)) .* ...
                g(delta, beta);
            
            % Define waypoints simple
%             waypoints = [ ...
%                 -Inf, Inf, ...
%                 ];
            % Choose waypoints intense.
            zz = z_bar(:)';
            waypoints = [ ...
                -Inf, Inf, ...
                beta(1), ...
                beta(1) - 10*beta(2), beta(1) + 10*beta(2), ...
                beta(1) - 10*min(sigma(:)), beta(1) + 10*min(sigma(:)), ...
                beta(1) - 10*max(sigma(:)), beta(1) + 10*max(sigma(:)), ...
                0, ...
                zz - 10*beta(2), zz + 10*beta(2), ...
                zz - 10*min(sigma(:)), zz + 10*min(sigma(:)), ...
                zz - 10*max(sigma(:)), zz + 10*max(sigma(:)), ...
                0, ...
                - 10*beta(2), 10*beta(2), ...
                - 10*min(sigma(:)), 10*min(sigma(:)), ...
                - 10*max(sigma(:)), 10*max(sigma(:)), ...
                ];            
            f = Twee.robust_integral(integrand, waypoints) - max(m, 0);
        end
    end
end
