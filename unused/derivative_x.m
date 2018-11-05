function  df_x = derivative_x(f,lambda,x)
% function  df_x = derivative_x(f,lambda,x)
%
% INPUT
% f       coefficients_polynomial
% lambda  scalar
% x       vector, length(x) == f.size_vector_field
%
% OUTPUT
% df_x    matrix f.size_vector_field*f.size_vector_field, derivative on x of f(lambda,x)


if length(x)~=f.size_vector_field
    error('dimensions do not match')
end

df_x = zeros(f.size_vector_field);
e = eye(f.size_vector_field);
for i = 1:f.size_vector_field
    for j =1:f.size_vector_field
        for n = 1: f.non_zero_el(i)
            if (f.powers{i}(n,j))>0
                df_x(i,j) = df_x(i,j) + ...
                    f.coefficients{i}(n) * lambda^f.powers_lambda{i}(n) ...
                    * (f.powers{i}(n,j)) * prod(horiz(x).^(f.powers{i}(n,:)-e(j,:)));  
            end
        end
    end
end
