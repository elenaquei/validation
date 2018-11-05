function [alpha_cont,s]=new_alpha_lorenz(s,nodes)
% OUTDATED
% function [alpha_cont,s]=new_alpha_lorenz(s,nodes)
%
% INPUT
% DIM      number of correlated Lorenz systems
% s        string of name of saved file
% nodes    number of nodes
% OUTPUT
% s        name of the saved system
%
% alpha for lorenz system
% 
% \dot x = sigma( y - z)
% \dot y = x(pho - z) - y
% \dot z = xy - beta z

if nargin<2 || isempty(s)
    s=sprintf('lorenz');
end

num_vec=3;
sigma=10;beta=8/3;


n_term =[3,4,3];
power_scal = cell(num_vec,1);
power_vector = cell(num_vec,1);
value = cell(num_vec,1);%{n_equations}(n_terms)
dot = cell(num_vec,1);
delay = dot;
for i = 1:num_vec
    dot{i} = cell(n_term(i),1);
    delay{i} = dot{i};
    for j = 1:n_term(i)
        dot{i}{j} = zeros(num_vec,1);
    end
    dot{i}{1}(i) = 1;
    power_vector{i} = cell(n_term(i),1);
    power_scal{i} = zeros(2,n_term(i));
    power_scal{i}(1,2:end)=1;
end
x = [1,0,0].';
y = [0,1,0].';
z = [0,0,1].';
i=0;
    power_scal{i*3+2}(2,2)=1;
    
    power_vector{i*3+1}{1} = zeros(num_vec,1);
    power_vector{i*3+1}{2} = zeros(num_vec,1);
    power_vector{i*3+1}{3} = zeros(num_vec,1);
    
    power_vector{i*3+2}{1} = zeros(num_vec,1);
    power_vector{i*3+2}{2} = zeros(num_vec,1);
    power_vector{i*3+2}{3} = zeros(num_vec,1);
    power_vector{i*3+2}{4} = zeros(num_vec,1);
    
    power_vector{i*3+3}{1} = zeros(num_vec,1);
    power_vector{i*3+3}{2} = zeros(num_vec,1);
    power_vector{i*3+3}{3} = zeros(num_vec,1);
    
    coord = i*3+(1:3);
    power_vector{i*3+1}{2}(coord) = y;
    power_vector{i*3+1}{3}(coord) = x;
    
    power_vector{i*3+2}{2}(coord) = x;
    power_vector{i*3+2}{3}(coord) = x+z;
    power_vector{i*3+2}{4}(coord) = y;
    
    power_vector{i*3+3}{2}(coord) = x+y;
    power_vector{i*3+3}{3}(coord) = z;
    
    value{i*3+1} = [1,-sigma,sigma];
    value{i*3+2} = [1,-1,1,1,0];
    value{i*3+3} = [1,-1,beta];


alpha_vec = polynomial_coefs(2, num_vec, num_vec, ...
    n_term, value,power_scal,power_vector, dot);

pol = polynomial_coefs(2,num_vec,0,zeros(0),cell(0),cell(0),cell(0),cell(0));
lin_coef = cell(3,1);
lin_coef{1} = ones(1,2);
lin_coef{2} = ones(1,num_vec,nodes*2+1);
lin_coef{3} = [0];
alpha_scal = scalar_eq(1,0,2,num_vec,lin_coef,pol);

alpha_cont = full_problem(alpha_scal, alpha_vec);

if nargout ==2
    save(s,'alpha_cont');
end