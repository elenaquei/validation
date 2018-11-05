function df_x = derivative_x_fourier(f,j)
% function df_x = derivative_x_fourier(f,j)
% 
% INPUT
% f       coefficients_polynomial
% lambda  scalar
% x       vector, length(x) == f.size_vector_field
% j       structure:
%     j.a       amplitude of periodic solution
%     j.omega   period of periodic solution 
%     j.y       Fourier series for y1 and y2
%     j.lambda  scalar parameter
%     j.x       vector stationary solution
%
% OUTPUT
% df_x    complex matrix, size 2 + f.size_vector_filed + (2n+1) * f.size_vector_field 

if length(j.x) ~= f.size_vector_field
    error('Dimensions do not agree');
end

n_nodes = (length(j.y) -1)/2;
M = f.size_vector_field;

[c, c_lambda]= compute_c_ik (f,j.x,j.lambda); % beeing c(k, iter) the coefficient d^i/dx^i f(x,lambda) and c_lamdba = d/dlambda c

%df_x = zeros( 2 + (2*n_nodes+1) * f.size_vector_field );

%df_x (3:3+(2*n_nodes+1)) =1;

J = repmat(vert(-n_nodes:n_nodes), f.size_vector_field,1);
Y = reshape( j.y.',1,[]).';

%df_x (3:end,1) = 1i * J .*Y;

e = eye(f.size_vector_field);

% da = zeros((2*n_nodes+1)*f.size_vector_field,1);

dx = zeros((2*n_nodes+1)*f.size_vector_field,M);
dy = zeros( (2*n_nodes+1)*f.size_vector_field);
dlambda = zeros( (2*n_nodes+1)*f.size_vector_field,1);
dlambda_small = zeros(f.size_vector_field,1);
dx_small = zeros(f.size_vector_field);

for k = 1: f.size_vector_field
    dyf_kn=0; % derivative of g with respect to y 
    dxf_kn = zeros((2*n_nodes+1),1); % derivative of gk with respect to xn
    dlambdaf_k = zeros((2*n_nodes+1),1); % derivative of gk with respect to lambda
    for n = 1: f.size_vector_field
        for i_iter = 1: prod(f.deg_f)
            
            i_index = convert(i_iter, f.deg_f);
            i_index_n = vert(i_index) + e(:,n);
            i_iter_n = convert_back(i_index_n, f.deg_f);
            
            dyf_kn = dyf_kn + i_index(n) * vert(c(k,i_iter) *  convolution(j.y,vert(i_index)-e(:,n)) * j.a^(sum(i_index)-1));
            if  (i_iter_n <=size(c,2))
            dxf_kn = dxf_kn + vert(c(k,i_iter_n) *  convolution(j.y,i_index) * j.a^(sum(i_index)-1));
            end
           
        end
       dy((k-1)*(2*n_nodes +1)+(1:(2*n_nodes+1)),(n-1)*(2*n_nodes +1)+(1:(2*n_nodes+1))) = toeplitz(dyf_kn);
       dx((k-1)*(2*n_nodes +1)+(1:(2*n_nodes+1)),n) = dxf_kn;
    end
    
    for i_iter = 1: prod(f.deg_f)
        i_index = convert(i_iter, f.deg_f);
        %da((n-1)*(2*n_nodes +1) +(1:(2*n_nodes+1))) = da((n-1)*(2*n_nodes +1) +(1:(2*n_nodes+1))) + ...
        %    sum(i_index) * vert(c(k,i_iter) *  convolution(j.y,i_index) * j.a^(sum(i_index)-2));
         dlambdaf_k = dlambdaf_k + vert(c_lambda(k,i_iter) *  convolution(j.y,i_index) * j.a^(sum(i_index)-1));
    end
    
    dlambda ((k-1)*(2*n_nodes +1)+(1:(2*n_nodes+1))) = dlambdaf_k;
    
    
    
    dlambda_k_small = 0; % derivative of tiny f with respect to lambda
    for i = 1: f.non_zero_el(k)
        if f.powers_lambda{k}(i)>0
            dlambda_k_small = dlambda_k_small + f.coefficients{k}(i)* (j.lambda).^(f.powers_lambda{k}(i)-1)*prod(horiz(j.x).^(f.powers{k}(i,:)));
        end
    end
    dlambda_small(k) = dlambda_k_small;
    
    for n = 1: f.size_vector_field
        dxf_kn_small = 0; % derivative of tiny f with respect to x
        for i = 1: f.non_zero_el(k)
            if all((f.powers{k}(i,:)-e(n,:))>=0)
            dxf_kn_small = dxf_kn_small + f.coefficients{k}(i)* f.powers{k}(i,n)*(j.lambda).^(f.powers_lambda{k}(i))*prod(horiz(j.x).^(f.powers{k}(i,:)-e(n,:)));
            end
        end
        dx_small(k,n) = dxf_kn_small;
    end
end

%J1 = 0*J;
%J1(1:2*n_nodes+1)=1;

p0_long_vec = reshape( repmat(horiz(f.phase_condition), 2*n_nodes+1,1), 1,[]);
p1_long_vec =reshape( repmat(horiz(f.amplitude_condition), 2*n_nodes+1,1), 1,[]);

zero_x = horiz(0*j.x);
%df_x = [0, 0, horiz(J1);
%    0, da(n_nodes+1), dy(n_nodes+1,:);
%    1i*vert(J.*Y), -da, 1i*j.omega*diag(J)-dy];


df_x = [ 0, 0, zero_x, p0_long_vec;
    0, 0, zero_x, p1_long_vec;
    dlambda_small, vert(zero_x), dx_small, zeros(length(j.x), length(J));
    -dlambda, 1i*vert(J.*Y), -dx, 1i*j.omega*diag(J)-dy];


return
end


