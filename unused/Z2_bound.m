function Z2vector=Z2_bound(DHm,A,xBar,alpha)
% function Z2vector=Z2_bound(DHm,A,xBar,alpha)
% 
% INPUT:
% DHm     complex matrix, can also be called Adagger, is the approximate 
%         derivative in the numerical solution xBar (given);
% A       complex matrix, numerical inverse of DHm;
% xBar    Xi_vector, numerical solution. Theoretically, there is no need
%         of having the numerical solution in this function, the goal is to have
%         the data if the problem inherent the solution (the size, the number
%         of nodes, ecc.);
% alpha  coefs, exact data of the problem.
%
% OUTPUT:
% Z0_vector   a standard vector of length = size_vector+size_scalar, the
%             values of this vector come from the Z1 bound in the associated pdf (we
%             refer to rigorous numerics for analytic solution od D.E..... for the
%             notation, mainly)
%
% In this function the coefficient of r.^2 in the radii polynomials is
% calculated.
% REMARK:
% A maximum radius need to be fixed in advance in order to
% simplify the calculations, it will be set as a global variable in the
% main, calling it RAD_MAX. If this is not done, the maximum radius is
% fixed by default to 10.^-3.
%
% It holds
%     Z2vector(i) = sum_{j=1}^{M+N} ||A_ij||  || [DDH(x+z)bc]_j ||
% and
%     || [DDH(x+z)bc]_j || <= sum_d sum_k sum_l | alpha_d,j d_k (d-e_k)_l |
%       ( |x|_c +r )^(d-e_k-e_l)

global RAD_MAX;
global Display
global use_intlab;

if isempty(RAD_MAX)
    RAD_MAX=1;
end
RAD_MAX=1;

A_ij= component_matrix_norm(A,xBar);

DDH=Function_second_derivative(xBar,alpha,RAD_MAX); 
% remark: the output WILL/SHOULD be a Xi_vector because only the convolution terms are
% left. BUT by now it's just a plain and normal vector. 

if use_intlab
Z2vector=sup(A_ij*DDH);
else
Z2vector=(A_ij*DDH);
end

if Display
    fprintf('Z2 computed, %d\n',Z2vector);
end