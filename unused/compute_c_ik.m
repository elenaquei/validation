function [c, dlambda_c]= compute_c_ik (f,x,lambda)
%function [c, dlambda_c]= compute_c_ik (f,x,lambda)

c = zeros(length(x), prod(f.deg_f+1)-1);
dlambda_c = c;
M = length(x);
for k = 1:M
    for i_iter = 1: prod(f.deg_f+1)-1
        i_index = convert(i_iter, f.deg_f);
        [c(k,i_iter), dlambda_c(k,i_iter)] = derivative (f,x,lambda,k,i_index);
    end
end
return
end


function [der, der_lambda] = derivative (f,x,lambda,k,i_index)
%function [der, der_lambda] = derivative (f,x,lambda,k,i_index)
% compute the derivative
% d^i/(dx^i) f_k(lambda,x)    with i multiindex, x vector
% and
% d/(d lambda) d^i/(dx^i) f_k(lambda,x)  

der =0;
der_lambda = 0;
for j = 1: f.non_zero_el(k)
    if all(f.powers{k}(j,:) - horiz(i_index(:))>=0)
        term =1;
        %vec_for = find((f.powers{k}(j,:) - horiz(i_index(:)))>=0);
        for n = 1:length(x) %vec_for 
            K = f.powers{k}(j,n) : -1: (f.powers{k}(j,n) - i_index(n)+1);
            new_prod = prod( K) * x (n)^(f.powers{k}(j,n) - i_index(n));
            term = term * new_prod;
        end
    else
        term =0;
    end
    
    if f.powers_lambda{k}(j)>0
        term_lambda = term * f.coefficients{k}(j)*(f.powers_lambda{k}(j))*lambda^(f.powers_lambda{k}(j)-1);
    else
        term_lambda =0;
    end
    
    term = term * f.coefficients{k}(j)*lambda^(f.powers_lambda{k}(j));
    
    der_lambda = der_lambda + term_lambda;
    der = der+ term;
end


end