function F = evalute_f_fourier(f,j)
% function f_x = evalute_f_fourier(f,j)
%
% INPUT
% f        coefficient_polyomial
% j        structure:
%     j.x       stationary solution
%     j.lambda  scalar parameter
%     j.a       amplitude of periodic solution = 0 
%     j.omega   period 1/b
%     j.y       Fourier series for y1 and y2
% 
% OUTPUT
% F        complex vector, F(j)    -- refer to pdf for further explanation

if f.size_vector_field ~=length(j.x)
    error('dimensions not compatible')
end

c= compute_c_ik (f,j.x,j.lambda);

n_nodes = (length(j.y) -1)/2;

f_x =0*j.y;
for k = 1: f.size_vector_field
    for i_iter = 1: prod(f.deg_f)
        i_index = convert(i_iter, f.deg_f);
        f_x(k,:) = f_x(k,:) + horiz(c(k,i_iter) *  convolution(j.y,i_index) * j.a^(sum(i_index)-1));
    end
end
J = repmat(vert(-n_nodes:n_nodes), f.size_vector_field,1);
Y = reshape( j.y.',1,[]).';

reshape_f_x = reshape( f_x.',1,[]).';

f_x_tot = 1i * j.omega * J .* Y - reshape_f_x;

% until here all okay

if size(j.y,1) ~= length(j.x)
    j.y=j.y.';
end
sum_Y = sum(j.y.');

% small_f_x = f(x,lambda)
small_f_x = 0*j.x;
for k =1: length(j.x)
    for i = 1: f.non_zero_el(k)
        small_f_x (k) = small_f_x(k) + f.coefficients{k}(i) * (j.lambda)^(f.powers_lambda{k}(i))*prod(horiz(j.x) .^ f.powers{k}(i,:));
    end
end


F = zeros(length(f_x_tot)+2+ f.size_vector_field,1);

F(1) = horiz(sum_Y)*vert(f.phase_condition);
F(2) = horiz(sum_Y)*vert(f.amplitude_condition) - 1 ;
F(2+(1: f.size_vector_field)) = vert(small_f_x);
F(3+ f.size_vector_field:end) = f_x_tot;


return
end
