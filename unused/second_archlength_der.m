function x_prime_prime = second_archlength_der(H,H0,H1,DxxHyy,x,x_prime,x_y_z)
% function x_prime_prime = second_archlength_der(H,H0,H1,DxxHyy,x,x_prime,x_y_z)
%
% compute the second derivative with respect to the archlength
%
% INPUT
% H         full_problem, at the endpoint we are considering
% H0, H1    full_problem, at the two endpoints of the segment
% DxxHyy    full_problem, DxxH*y*y, input dimentsions 3*input dim of H0
% x         Xi_vector, fitting input for H0, H1
% x_prime   Xi_vector, fitting input for H0, H1, first derivative
% x_y_z     Xi_vector, fitting input for DxxHyy
%
% OUTPUT
% x_prime_prime     vector, second derivative with respect to the
%                   archlength parameter

DxH = derivative_to_matrix(derivative(H,x,0));
v_0 = extract_all_lin_coef(H0.scalar_equations);
v_1 = extract_all_lin_coef(H1.scalar_equations);
v_Delta = v_1 - v_0;


DsxH = v_Delta;
DsxH = cat(1,DsxH, zeros(length(v_0)-H.scalar_equations.number_equations_lin,length(v_0)));
DxxH_yy = Xi_vec2vec(apply(DxxHyy,x_y_z,0));

x_prime_prime = (DxH \ (-2*DsxH* Xi_vec2vec(x_prime)-DxxH_yy));
