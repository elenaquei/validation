function [flag, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F0,F1, x0,x1,num_variable)
% function [flag, x_prime0, x_prime1,x_prime_prime0,x_prime_prime1] = if_saddle_numerical(F0,F1, x0,x1,num_variable)
%
% numerical test for a saddle node in the segment [x0,x1].
%
% INPUTS
% F0, F1        full_problem, at the beginning and at the end of the
%               segment
% x0, x1        Xi_vector, numerical solution
%num_variable   integer, position of the variable to check (DEFAULT: all)
%
% OUTPUTS
% flag          bool, 1 if a saddle was found numerically
% x_prime0, x_prime1  vectors, approximate first derivatives w.r.t. the
%               archlength
% x_prime_prime0, x_prime_prime1   vectors, (very) approximate second
% derivativec w.r.t. archlength

if nargin<5 || ~exist('num_variable','var')
    num_variable = 1:x0.size_scalar;
end


flag = 0;

Dx0H = derivative_to_matrix(derivative(F0,x0,0));
Dx1H = derivative_to_matrix(derivative(F1,x1,0));

x0_vec = vert(Xi_vec2vec(x0));
x1_vec = vert(Xi_vec2vec(x1));
x_Delta = x1_vec - x0_vec;
v_0 = extract_all_lin_coef(F0.scalar_equations);
v_1 = extract_all_lin_coef(F1.scalar_equations);
v_Delta = v_1 - v_0;

c_0 = vert(F0.scalar_equations.linear_coef{3});
c_1 = vert(F1.scalar_equations.linear_coef{3});
c_Delta = c_1-c_0;

Ds0H = v_Delta*x0_vec + 0*v_0*x_Delta +c_Delta;
Ds0H = cat(1,Ds0H, zeros(length(v_0)-F0.scalar_equations.number_equations_lin,1));
Ds1H = v_Delta*x1_vec + 0*v_1*x_Delta +c_Delta;
Ds1H = cat(1,Ds1H, zeros(length(v_0)-F0.scalar_equations.number_equations_lin,1));

DsxH = v_Delta;
DsxH = cat(1,DsxH, zeros(length(v_0)-F0.scalar_equations.number_equations_lin,length(v_0)));
% it was like this, why? 
% Ds0H = dot(v_Delta,x0_vec) + dot(x_Delta, v_0) +c_Delta;
% Ds1H = dot(v_Delta,x1_vec) + dot(x_Delta, v_1) +c_Delta;
% it doesn't make sense, why?

x_prime0 = - Dx0H \ Ds0H;
x_prime1 = - Dx1H \ Ds1H;
x_prime_prime0 = -2* (Dx0H \ (DsxH* x_prime0));
x_prime_prime1 = -2* (Dx1H \ (DsxH* x_prime1)); 

% x_prime0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_prime0,x0)));
% x_prime1 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_prime1,x0)));
% x_prime_prime0 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_prime_prime0,x0)));
% x_prime_prime1 = Xi_vec2vec(symmetrise(vec2Xi_vec(x_prime_prime1,x0)));

if any(real(x_prime0(num_variable)).*real(x_prime1(num_variable))<0)
%     if all(abs(real(x_prime0(num_variable)).*real(x_prime1(num_variable)))<10^-10)
%         warning('Impossible to expect validation, boundaries too close to zero')
%         flag=0;
%     else
        flag = 1;
%     end
end

return
end

% function v = extract_coefficients(F)
% % function v = extract_coefficients(F)
% % 
% % takes F and extract the coefficients of the linear scalar equations
% scal_eq = F.scalar_eq;
% v = zeros(scal_eq.size_scalar + scal_eq.size_vector *(2*scal_eq.nodes +1),scal_eq.number_lin_eq); % check naming
% v(1:scal_eq.size_scalar,:) = scal_eq.lin_coef{1};
% v(scal_eq.size_scalar+1:end,:) = reshape(scal_eq.lin_coef{2}(:,:,:),[],scal_eq.number_lin_eq) ; % check the ordering!
% end