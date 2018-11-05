function f_x = evalute_f(f,lambda,x)
% function f_x = evalute_f(f,lambda,x)
%
% INPUT
% f        coefficient_polyomial
% lambda   scalar
% x        vector, length=f.size_vector_field
% 
% OUTPUT
% f_x      f(lambda,x)

if f.size_vector_field ~=length(x)
    error('dimensions not compatible')
end

f_x =0*x;
for k = 1: f.size_vector_field
    for  n = 1:f.non_zero_el(k)
        if any (f.powers{k}(n,:)>0)
            f_x(k) = f_x(k) + f.coefficients{k}(n)*lambda^(f.powers_lambda{k}(n))*prod(horiz(x).^f.powers{k}(n,:));
        else
            f_x(k) = f_x(k) + f.coefficients{k}(n)*lambda^(f.powers_lambda{k}(n));
        end
    end
end


