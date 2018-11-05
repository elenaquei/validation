function [alpha_coef, coef, x_Xi] = converting_data(j,f)
% function [alpha_coef, coef, x_Xi] = converting_data(j,f)
%
% INPUT
% j     structure:
%     j.x       stationary solution
%     j.lambda  Hopf bifurcation parameter
%     j.a       amplitude of periodic solution = 0 
%     j.omega   period 1/b
%     j.y       Fourier series for y1 and y2
% f         coefficients_plynomial
% 
% OUTPUT - compatible with old code
% alpha      of type 
% coef       cell{3} to store G, linear - REMEMBER TO OVERWRITE THE CALLS
%            TO G
% x_Xi       Xi_vector



% alpha_cont = coefs(); % overwritten later on 
coef = cell(3,1);
%x_Xi = Xi_vector();
size_vec = length(x);

% Xi_vec = vec2Xi_vec(vec,size_scal,size_vec,nodes);
% alpha=coefs(Sscal,Svec,Dscal,Dvec,Non0,Pscal,Pvec,val);
% coef{1} = [p0, p1, 0] (0 to be overwritten by f(x,lambda)

%% x_Xi start
% check y vertical
if size(j.y,1) == length(x)
    j.y=j.y.';
elseif size(j.y,2) ~= length(x)
    error('x and y need to be compatible');
end

% turn y in a vector
y_vec = reshape(j.y,1,numel(j.y));

vec = [j.omega, j.lambda, j.a, horiz(j.x), y_vec];

nodes = (size(j.y,1)-1)/2;

if ~isinteger(nodes)
    error('number of nodes not integer');
end

x_Xi = vec2Xi_vec(vec,size_vec+3 ,size_vec,nodes);
% x_Xi completed

%% start coef
coef{3} = [0,1,zeros(1,size_vec)];
coef{2} = zeros(size_vec +2, size_vec, 2*nodes +1);
for i =1:size_vec
    coef{2}(1,i,:) = f.phase_condition(i);
    coef{2}(2,i,:) = f.amplitude_condition(i);
end
coef{1} = zeros(size_vec +2, size_vec+2);
% coef completed

%% start alpha

% HERE somethign smart needed! 
size_scal = size_vec + 3; 
% scal  = [ omega, lambda, x, a]
size_y = size_vec;
Non0= zeros(size_y,1);
Pscal = cell(size_y,1);
Pvec = cell(size_y,1);
val = cell(size_y,1);


Non0_non_cont= zeros(size_y,1);
Pscal_non_cont = cell(size_y,1);
Pvec_non_cont= cell(size_y,1);
val_non_cont = cell(size_y,1);
for j = 1 : size_y % equation we are considering
    counter = 0;
    Non0 = prod(prod(f.powers_vector{j}+1));
    Pscal{j} = zeros(Non0, size_scal);
    Pvec{j} = zeros(Non0, size_y);
    val{j} = zeros(Non0, 1);
    
    counter_non_cont =0;
    
    for m=1:f.non_zero_el(j)
        p = f.powers_vector{j}(m,:);
        for i = 2: prod(p+1)
            counter = counter+1;
            k = convert(i,p);
            K = construct_product_matrix(size_y,k,p);
            c = f.coefficients{j}(m);
            c_new= c*prod(prod(K));
            q_new = [f.powers_lambda{j}(m), horiz(p-k), sum(k)-1];
            p_new = horiz(k);
            
            Pscal{j}(counter,:) = q_new;
            Pvec{j}(counter,:) = p_new;
            val{j}(conter) = c_new;
            if sum(k)-1==0
                counter_non_cont = counter_non_cont +1;
                Pscal_non_cont{j}(counter_non_cont,:) = q_new(1:end-1);
                Pvec_non_cont{j}(counter_non_cont,:) = p_new;
                val_non_cont{j}(counter_non_cont) = c_new;
            end
        end
    end
    
    Non0_non_cont(j) = counter_non_cont;
end

alpha_coef=coefs_scal_cont(size_scal,size_y,Dscal,Dvec,Non0,Pscal,Pvec,val);

alpha_coef = compact_alpha(alpha_coef); % avoid doubling of powers


s = 'Hopf';
save(sprintf('%s_cont',s),'alpha_coef');

% alpha_cont completed

%% alpha_non_cont

alpha_coef=coefs_scalar(size_scal-1,size_y,Dscal,Dvec,Non0,Pscal_non_cont,Pvec_non_cont,val_non_cont);
save(s,'alpha_coef');


end




%% help function 
function K = construct_product_matrix(N, k,p)

K = ones(N, max(k));
for i = 1:max(k)
    index = find( (p-i+1).*(k-i)>0);
    K(index,i) = p(index)-i+1;
end
end


