function [x, error_out, residual] = Newton_scalar( f, lambda, x, max_iter, min_res)
% function x_new = Newton_scalar( f, lambda, x)
% 
% INPUT
% f         coefficients_polynomial
% lambda    fixed scalar
% x         firs approximation of the solution
% max_iter  maximum number of iterations (DEFAULT 20)
% min_res   requested residual (DEFAULT 10^-8)
%
% OUTPUT
% x_new     final approximation
% error     0 if failed, 1 if good, 2 if no need of iterations
% residual  final residual


if nargin<3
    error('not enough inouts');
end

default_iter = 20;
default_res = 10^-8;

if nargin<4 || isempty(max_iter)
    max_iter= default_iter;
end
if nargin<5 || isempty(min_res)
    min_res = default_res;
end
error_out =0;
f_x = evalute_f(f,lambda,x);
iter =0;
while norm(f_x)>min_res && iter<max_iter
    iter=iter+1;
    df_x = derivative_x(f,lambda,x);
    x = vert(x) - df_x \ vert(f_x);
    f_x = evalute_f(f,lambda,x);
    
    res(iter) = norm(f_x);
    
end
residual = norm(f_x);

% checking residual
if residual > min_res
    warning('Newton did not converge')
    return
end

% outputing 'ok' flag
if iter ==0
    error_out =2;
else
    error_out=1;
end
