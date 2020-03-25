% This script creates the following three Matlab function files:
%
% 1)  t_distribution.m  : t_distribution(x,[m,s,v]) evaluates the p.d.f 
%                         of the random variable X = m + st_v where t_v is
%                         t-distributed with v degrees of freedom.
%
% 2)  d_t_distribution.m: d_t_distribution(x,[m,s,v]) reports the gradient 
%                         of the p.d.f of the random variable X = m + st_v 
%
% 3) dd_t_distribution.m: dd_t_distribution(x,[m,s,v]) reports the hessian 
%                         of the p.d.f of the random variable X = m + st_v
% 
% NOTE: The gradient and derivatives are generated using the Symbolic Math 
%       toolbox.
%       The functions are saved in the main directory under
%       Matlab/Functions
%
% This version: November 29th, 2017.

clear;
reset(symengine);

%% 1) Construct the symbolic variables
syms x m s v;

assume(x, 'real')

assume(m, 'real');

assume(s > 0);

assume(v > 0);

%% 2) Construct the symbolic functions

t_distribution = ...
                1 / s ...
                * gamma((v+1)/2) ...
                / sqrt(v*pi) / gamma(v/2) ...
                * (1 + ((x / s - m / s)^2) / v)^(-(1+v)/2);

t_gradient     = gradient(t_distribution, [m, s, v]);

t_hessian      = hessian(t_distribution, [m, s, v]);

%% 3) Construct the Matlab function files
matlabFunction(t_distribution, ...
  'Vars', {x, [m, s, v]}, ...
  'file', './matlab/functions/t_distribution.m');

matlabFunction(t_gradient, ...
  'Vars', {x, [m, s, v]}, ...
  'file', './matlab/functions/d_t_distribution.m');

matlabFunction(t_hessian, ...
  'Vars', {x, [m, s, v]}, ...
  'file', './matlab/functions/dd_t_distribution.m');

% clear;
% reset(symengine);