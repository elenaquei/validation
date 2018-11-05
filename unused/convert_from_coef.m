function alpha = convert_from_coef(alpha_coef)
% function alpha = convert_from_coef(alpha_coef)
%
% INPUT
% alpha_coef     instance of the class coefs
% OUTPUT
% alpha          instance of the class polynomial_coefs
%
% alpha is equivalent to alpha_coef
alpha = polynomial_coef();

alpha.size_scalar=alpha_coef.size_scalar;
alpha.size_vector= alpha_coef.size_vector;
alpha.deg_scalar = alpha_coef.deg_scalar;
alpha.deg_vector = alpha_coef.deg_vector;
alpha.n_equation = alpha_coef.size_vector;
alpha.n_terms = alpha_coef.non_zero_el;
alpha.delay = get_default_delay(alpha);
alpha.monomial_per_term = get_default_monomial(alpha);
alpha.dot = get_default_derivative(alpha);
alpha.value = alpha_coef.value;
end