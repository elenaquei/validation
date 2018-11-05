% function test_case
% test problem %
%
% in this script, functions are tested one by one, on a symple problem,
% such as 
%    x' = - lambda y
%    y' =   lambda x
%    x0 = 1
%    y0 = 0
% with solution 
%    x = cos ( lambda t )
%    y = sin ( lambda t )

global use_intlab
use_intlab=0;

m=20;

value=cell(2);
value{1}=[1i, 1/(2*pi)];
value{2}=[1i, -1/(2*pi)];
powers_scalar=cell(2);%zeros(2,2,2);
powers_scalar{1}=[0,1;0,1];
powers_scalar{2}=[0,1;0,1];

powers_vector=cell(2);
powers_vector{1}=[1,0;0,1];
powers_vector{2}=[0,1;1,0];
non_zero=cell(2);
non_zero{1}=2;non_zero{2}=2;
alpha_coef=coefs(2,2,2,1,non_zero,powers_scalar,powers_vector,value);
save('test_case_cont','alpha_coef');


value=cell(2);
value{1}=[1i, 1/(2*pi)];
value{2}=[1i, -1/(2*pi)];
powers_scalar=cell(2);
powers_scalar{1}=[0,1];
powers_scalar{2}=[0,1];

powers_vector=cell(2);
powers_vector{1}=[1,0;0,1];
powers_vector{2}=[0,1;1,0];
alpha_coef=coefs(1,2,2,1,non_zero,powers_scalar,powers_vector,value);

save('test_case','alpha_coef');

Name_system='test_case';
flag_plot=0;
n_nodes=20;
maxiter=100;
min_res=10^-14;
delta=0.1;

[x0, x_dot0, ~,xShort0,alpha_coef,coefs_linear] = ...
    compute_solution_and_derivative(Name_system,flag_plot,n_nodes,maxiter,min_res);

% running until here, now would be cool to keep on going! 

use_intlab=0;

[x1,x_dot1,DH1,coefs_linear1] = update(x0,delta,x_dot0,...
    alpha_coef,coefs_linear,maxiter,min_res);

use_intlab = temp_intlab;


