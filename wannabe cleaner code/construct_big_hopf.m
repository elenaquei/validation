function big_Hopf = construct_big_hopf( ...
       lambda0, x0, f, n_nodes, eigenvec, eigenval)
% function [s, x_n] = continuation_Hopf ( lambda0, x0, f, ....
%       n_nodes, n_iter, eigenvec, eigenval, sign_FLC)
%
% INPUTS
% lambda0       paramenter of the Hopf bifurcation
% x0            scalar values of the Hopf bifurcation
% f             vector field defining the problem
% n_nodes       integer, number of nodes to use for Fourier
% n_iter        number of validated continuation iterations to do (DEFAULT:
%               100)
% h             step size in the continuation code
% s             string, path where to save results
% eigenvec      approximation of the eigenvector at the Hopf bifurcation
% eigenval      approximation of the eigenvalue at the Hopf bifurcation
%
% OUTPUTS
% s              string, path of saved solutions
% x_n            Xi_vector, last validated solution


if length(x0)<2
    error('Need at least 2 dimensions to have a Hopf bifurcation')
end

% transformation of the vector field into the bigger system needed for Hopf
big_Hopf = Taylor_series_Hopf(f,n_nodes);


% analytic contruction of the solution at the Hopf bifurcation
T_star = imag(eigenval);
a_star = 10^-4;
if nargin>9
    lambda0 = lambda0+ sign_FLC*a_star^2;
end

% guess of an approximation of the solution (assuming the eigenvalues and
% eigenvectors at the Hopf bifurcation to give a good enough idea 
y = zeros(length(x0),2*n_nodes+1);

y (:,n_nodes)= conj(eigenvec);
y (:,n_nodes+2)= eigenvec;

% constructing the full solution of the blowed up system
sol = Xi_vector([1/T_star, lambda0, a_star, horiz(x0)],y);

% adding 2 scalar condition
% first scalar condition: integral condition
big_Hopf.scalar_equations =fancy_scalar_condition(sol,big_Hopf.scalar_equations,1);

% second sclar condition: some constraint on the norm of z (the rescaled
% periodic orbit)
p1 = Xi_vec2vec(sol).';
p1(1:sol.size_scalar) = 0;
p1 = conj(p1);
const_p1 = -sum(p1.*conj(p1));
big_Hopf.scalar_equations = change_lin_coef_vector(big_Hopf.scalar_equations,[p1.';const_p1],2);

return
end