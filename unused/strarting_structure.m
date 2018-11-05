function [j,f] = strarting_structure(f, x, lambda, sin_cos, b, n_nodes)
%function [j,f] = strarting_structure(f, x, lambda, sin_cos, b, n_nodes)
%
% INPUT
% f         coefficients_plynomial
% x         stationary solution
% lambda    hopf bifurcation coefficient
% sin_cos   boolean if y is (sint,cost) or (cost, sint)
% b         imaginary part of the eigenvalue of f in the Hopf bifurcation
% n_nodes   number of nodes of the output
%
% OUTPUT
% j     structure:
%     j.x       stationary solution
%     j.lambda  Hopf bifurcation parameter
%     j.a       amplitude of periodic solution = 0 
%     j.omega   period 1/b
%     j.y       Fourier series for y1 and y2
% f      coefficients_polynomial, now added
%     p0    vector, for the first scalar equation
%     p1    vector, for the second scalar equation

sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = -1i/2; sin_four(n_nodes+2) = 1i/2;
if sin_cos
    j.y = [sin_four;cos_four];
    p0 = [1,0];
    p1 = [sum(sin_four), sum(cos_four)];
else 
    j.y = [cos_four;sin_four];
    p0=[0,1];
    p1 = [sum(cos_four), sum(sin_four)];
end
j.a = 0;
j.omega = 1/b;
j.x = x;
j.lambda = lambda;
f = scalar_condition(f, p0, p1);

return 
end