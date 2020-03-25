addpath('./functions');

a(1) = tpdf_with_gradient(3, [0 1 7]);
b(1) = tpdf(3, 7);

m = 5;
s = 3;
v = 2;
x = -1;
a(2) = tpdf_with_gradient(x, [m s v]);
b(2) = tpdf((x-m)/s, v) / s;

if norm(a-b) > 10^(-10)
    error('t distribution problem');
end

clear;
[g, dg] = tpdf_with_gradient([-3, 1], [1 2 3]);

if ~isequal(size(g), [1, 2]) || ~isequal(size(dg), [3, 2])
    error('gradient has the wrong size!');
end

disp('test_t_distribution passed');
