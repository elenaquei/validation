% validation of saddle nodes in a standard ODE system 

% the validation will follow the outlined path:
% 1) find x0 and x1 
% 2) extend the system around them 
% 3) compute the extended system - not at all straightforward Newton
% 4) validate the extended system and solution - automatic with the code as
% it is now

%% find the saddle node
%
% given an ODE +scalar eq.s system 
% H(l,v)=[g(l,x)  
%         \dot x - f(l,x)] =0 
%
% 1-D underdetermined 

% EXAMPLE: Rychkov
%
% x_dot = y - x^5 + x^3 - lambda x
% y_dot = -x

global first_run
global nu
global use_intlab 
%global talkative 
global RAD_MAX
%talkative = 5;
use_intlab = 0;
nu = 1.01;
RAD_MAX = 10^-2;

if isempty(first_run)
    addpath(genpath('../'))
    addpath(genpath('../../'))
    startintlab
    first_run =1;
end

n_nodes = 80;

% load('fifth_order_validation_10','x0')
% lambda = x0.scalar(2);
% x0_load = x0;
% x0_load.size_scalar = 1;
% x0_load.scalar = [x0.scalar(1)/(2*pi)];
% 
% sin_four = zeros(1,2*n_nodes+1);
% cos_four = sin_four;
% cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
% sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
% y_vec = 0.8116 * [sin_four; cos_four];
% sol = Xi_vector([1, lambda],y_vec);
% sol_fixed = Xi_vector(1,y_vec);
% 
% 
% vector_field_string = '-dot x1 + l1 x2 - l1 x1^5 + l1 x1^3 - l1l2 x1 \n -dot x2 - l1 x1';
% vector_field_string_fixed_l = '-dot x1 + l1 x2 - l1 x1^5 + l1 x1^3 - lambda l1 x1 \n -dot x2 - l1 x1';
% vector_field_string_fixed_l = strrep(vector_field_string_fixed_l, 'lambda' , num2str(lambda,16)); 
% 
% f_1 = @(x,y) y - conv(x,conv(x,conv(x,conv(x,x,'same'),'same'),'same'),'same') + conv(x,conv(x,x,'same'),'same') - lambda *x;
% f_2 = @(x,y) - x;
% F = @(X) [f_1(X(1,:),X(2,:)); f_2(X(1,:),X(2,:))];
% 
% 
% %function x = Newton_silly(x0)
%     x = x0_load.vector;
%     x = [1,x(1,:),x(2,:)].';
%     middle = @(x) length(x)/2;
%     split = @(x)[horiz(x(1:middle(x)));horiz(x(middle(x)+1:end))];
%     first_in = @(x) x(1,:);
%     first = @(y) first_in(split(y));
%     F_t = @(t,y) [sum(first(y(2:end))); y(1)*reshape(F(split(y(2:end))),[],1)];
%     DF = @(x) numjac(F_t,0,x, F_t(0,x), 10^-5);
%     F_newton = @(x) F_t(0,x);
%     for i = 1:10
%         x = x - DF(x)\F_newton(x); 
%     end
% %end
% 
% 
% 
% 
% 
% 
% f_ode = @(t,x) [ x(2) - x(1)^5 + x(1)^3 - lambda * x(1) ; - x(1)];
% 
% [~,y] = ode45(f_ode,[0:0.01:40*pi], [0.8,0]);
% [t,y] = ode45(f_ode,[0:0.01:2*pi], y(end,:));
% x0_fixed = time_series2Xi_vec(t,y,n_nodes);
% x0 = x0_fixed;
% x0.size_scalar = 2;
% x0.scalar = [x0.scalar, lambda];
% 
% f_fixed = from_string_to_polynomial_coef(vector_field_string_fixed_l);
% f = from_string_to_polynomial_coef(vector_field_string);
% 
% scalar_eqs_fixed = default_scalar_eq(x0_fixed,1);
% F_fixed = full_problem(scalar_eqs_fixed,f_fixed);
% 
% F_fixed.scalar_equations.linear_coef{2}(1,1,:)=0;
% F_fixed.scalar_equations.linear_coef{2}(1,2,:)=1;
% F_fixed.scalar_equations.linear_coef{3} =0;
% 
% x0_N = Newton_2(x0_load,F_fixed,[],0.009); % still very bad: period 10^6
% 
% x_dot = kernel(derivative_to_matrix(derivative(F,x0_rand,0)));
% const = -Xi_vec2vec(x0_rand).'*x_dot;
% F.scalar_equations = change_lin_coef_vector(F.scalar_equations,[x_dot;const],2); 
% x0_N = Newton_2(x0_rand,F);
% 
% scalar_eqs = default_scalar_eq(sol,1);
% scalar_eqs.linear_coef{2}(1,1,:)=0;
% scalar_eqs.linear_coef{2}(1,2,:)=1;
% scalar_eqs.linear_coef{3} = 0.6;
% F = full_problem(scalar_eqs, f);
% continuation(sol,F)

n_iter = 10;
h = 5*10^-4;
s = 'saddle_test';
s_num = 'saddle_num';
min_res_N = 10^-9;
mu = 2.5;
f = @(t,x) [ (x(2)-(x(1).^5 - mu *x(1)^3+x(1))); -x(1)];
[t,y]=ode45(f,[0:0.01:6.8],[1.3;-0.8]);
x0 = time_series2Xi_vec(t,y,n_nodes);

string_vf = 'dot x1 - l1  x2 + l1 x1^5 - 2.5 l1 x1^3 + l1 x1 \n dot x2 + l1 x1';
f = from_string_to_polynomial_coef(string_vf);

scalar_eqs_fixed = fancy_scalar_condition(x0);
% scalar_eqs_fixed = default_scalar_eq(x0,1);
% scalar_eqs_fixed.linear_coef{3} = -1;
F_fixed = full_problem(scalar_eqs_fixed,f);
x0_N = Newton_2(x0,F_fixed,[],min_res_N);

string_vf = 'dot x1 - l1  x2 + l1 x1^5 - l2 l1 x1^3 + l1 x1 \n dot x2 + l1 x1';
f = from_string_to_polynomial_coef(string_vf);
x0_N.size_scalar = 2;
x0_N.scalar = [x0_N.scalar,mu];

scalarEqs =  fancy_scalar_condition(x0_N); %default_scalar_eq(x0_N,1);
F_not_square = full_problem(scalarEqs,f);

DF_not_square = derivative_to_matrix(derivative(F_not_square,x0_N,0));
x_dot = kernel(DF_not_square);
[~,index] = max(abs(real(x_dot(1:x0_N.size_scalar))));
angle = atan( imag(x_dot(index))/real(x_dot(index)));
x_dot = exp( - 1i * angle) * x_dot;
x_dot = Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot,x0_N)));
if x_dot(2)>0
    x_dot = -x_dot;
end

h_num = 0.5*10^-3;
[s_num,x1] = continuation_numerical( x0_N, F_not_square, 135, h_num, x_dot,'test_get_to_saddle', min_res_N, 0 ,1); %
load(s_num)
no_fancy_condition = 0;
no_Hopf_system = 0;
check_for_saddle =1;
[s, x_n] = continuation ( x1, F_not_square, 22, h, x_dot_1,'test_validate_saddle', min_res_N,no_Hopf_system,no_fancy_condition,check_for_saddle);
