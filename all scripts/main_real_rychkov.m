% validation of saddle nodes in a standard ODE system 
% EXAMPLE: Rychkov
%
% x_dot = y - x^5 + x^3 - lambda x
% y_dot = -x

% working as of 4th July 2018 (RESCALING IN bigger_saddle_system BY HAND)
% 300 intervals, rescaling 500


global nu
global use_intlab 
%global talkative 
global RAD_MAX
global rescaling_saddle
global refinement_saddle
rescaling_saddle = 11.5*10^5;% still too big
refinement_saddle =750;% from 500, good enough to find the saddle numerically
%talkative = 2;
use_intlab = 0;
nu = 1.01;
RAD_MAX = 10^-2;
% BIG IDEA: split results of radii polynomial approach into a vector, one
% bound for each element of x


try 
    intval(1);
catch
    addpath(genpath('../'));
    addpath(genpath('../../'))
    startintlab;
end

n_nodes = 44;

n_iter = 10;
h = 5*10^-4;
s = 'rychkov_test';
s_num = 'rychokov_num';
s_temp = 'temp';
min_res_N = 10^-9;
%mu = 2.5;
mu = 0.224;
%f = @(t,x) [ (x(2)-(x(1).^5 - mu *x(1)^3+x(1))); -x(1)];
f = @(t,x) [ (x(2)-(x(1).^5 - x(1)^3+mu *x(1))); -x(1)];

[t,y]=ode45(f,[0:0.01:6.8],[1.3;-0.8]);
[t,y]=ode45(f,[0:0.01:6.8],y(end,:));
x0 = time_series2Xi_vec(t,y,n_nodes);

%string_vf = 'dot x1 - l1  x2 + l1 x1^5 - 2.5 l1 x1^3 + l1 x1 \n dot x2 + l1 x1';
string_vf = 'dot x1 - l1  x2 + l1 x1^5 - l1 x1^3 + 0.224 l1 x1 \n dot x2 + l1 x1';
f = from_string_to_polynomial_coef(string_vf);

scalar_eqs_fixed = fancy_scalar_condition(x0);
% scalar_eqs_fixed = default_scalar_eq(x0,1);
% scalar_eqs_fixed.linear_coef{3} = -1;
F_fixed = full_problem(scalar_eqs_fixed,f);
x0_N = Newton_2(x0,F_fixed,[],min_res_N);

%string_vf = 'dot x1 - l1  x2 + l1 x1^5 - l2 l1 x1^3 + l1 x1 \n dot x2 + l1 x1';
string_vf = 'dot x1 - l1  x2 + l1 x1^5 - l1 x1^3 + l1 l2 x1 \n dot x2 + l1 x1';
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

h_num = 5*10^-3;
[s_num_temp,x1] = continuation_numerical( x0_N, F_not_square, 30, h_num, x_dot,s_temp, min_res_N, 0 ,1); %
load(s_num_temp)
[s_num,x1] = continuation_numerical( x1, F_not_square, 46, h_num, -x_dot_1,s_num, min_res_N, 0 ,1); %
load(s_num)
no_fancy_condition = 0;
no_Hopf_system = 0;
check_for_saddle =1;
[s, x_n] = continuation ( x1, F_not_square, 40, h, x_dot_1,s, min_res_N,no_Hopf_system,no_fancy_condition,check_for_saddle);% saddle at iteration 28
