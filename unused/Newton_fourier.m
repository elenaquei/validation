function [j, error_out, residual] = Newton_fourier ( f, j, max_iter, min_res)
% function [j, error_out, residual] = Newton_fourier ( f,  j, max_iter, min_res)
% 
% INPUT
% f          coefficients_polynomial
% j          initial guess, structure:
%     j.x       stationary solution
%     j.lambda  Hopf bifurcation parameter
%     j.a       amplitude of periodic solution
%     j.omega   period of periodic solution 
%     j.y       Fourier series for y1 and y2
% max_iter  maximum number of iterations (DEFAULT 20)
% min_res   requested residual (DEFAULT 10^-8)
%
% OUTPUT
% j          final approximation, structure:
%     j.x       stationary solution
%     j.lambda  Hopf bifurcation parameter
%     j.a       amplitude of periodic solution
%     j.omega   period of periodic solution 
%     j.y       Fourier series for y1 and y2
% error     0 if failed, 1 if good, 2 if no need of iterations
% residual  final residual
%
if nargin<2
    error('not enough inouts');
end

default_iter = 20;
default_res = 10^-8;

if nargin<3 || isempty(max_iter)
    max_iter= default_iter;
end
if nargin<4 || isempty(min_res)
    min_res = default_res;
end

error_out =0;
f_x = evalute_f_fourier(f,j);
iter =0;
while norm(f_x)>min_res && iter<max_iter
    iter=iter+1;
    df_x = derivative_x_fourier(f,j);
    vec_j = reshape_j(j);
    new_vec_j = vert(vec_j) - df_x \ vert(f_x);
    j = restructure(new_vec_j,j);
    f_x = evalute_f_fourier(f,j);
    
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

return
end

function vec_j = reshape_j(j)
% function vec_j = reshape_j(j)
%
% turns j into a vector (lambda, omega, a, x, y)
Y = reshape( j.y.',1,[]).';
vec_j = [j.lambda, j.omega,horiz(j.x), horiz(Y)].';


end


function j = restructure( vec_j, old_j)
%function j = restructure( vec_j, old_j)
j.a = old_j.a;
j.lambda = vec_j(1);
j.omega = vec_j(2);
M=length(old_j.x);
j.x = vec_j(2+(1:M));

j.y = reshape(vec_j(3+M:end),M,[]);

end
