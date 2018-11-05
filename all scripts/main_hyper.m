% fancy 4D cahotic
%
% still not working properly
%
% 
%

global first_run
global nu
global use_intlab 
%global talkative 
global RAD_MAX
%talkative = 1;
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;

if isempty(first_run)
    addpath(genpath('../'))
    startintlab
    first_run =1;
end
%clearvars

% b=c<sqrt(2e)
% d neq 0
% e>0
e = 2;
b=1;
c=b;
d=10;

% alpha_null=0;

string_hyper = '- dot x1 +a x1 -ax2-x2 x3+x4  \n - dot x2 -b x2 +x1 x3 \n - dot x3 -cx3+Dx1+x1x2 \n - dot x4 -ex1-ex2'; % fancy hyper
string_hyper_vars = strrep(string_hyper, 'b' , num2str(b)); % plugging in b
string_hyper_vars = strrep(string_hyper_vars, 'c' , num2str(c)); % plugging in c
string_hyper_vars = strrep(string_hyper_vars, 'D' , num2str(d)); % plugging in d
string_hyper_vars = strrep(string_hyper_vars, 'e' , num2str(e)); % plugging in e

string_hyper_cont = strrep(string_hyper_vars, 'a', 'l2'); % setting a as the second scalar variable


% some elements useful for the computation and the validation
n_nodes = 5; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 5; % number of iterations
step_size = 10^-3; % initial step size (then adapted along the validation
s = 'Hopf_hyper'; % where the solutions are stored

vectorfield = strrep(string_hyper_cont, 'l1' , '');
vectorfield = strrep(vectorfield, 'l2' , 'l1');
% string defining the vector field of the Hopf normal form 

f = from_string_to_polynomial_coef(vectorfield); % trasnformation into a vectorfield that can be used
fn = @fn_hyper;
% definition of the solution
%x_star 
%lambda_star 
% dummy = load('hopf_in_hyper');
% 
% % validated values converted to doubles 
% eigenval = conj(mid(dummy.eigenval));
% eigenvec = conj(mid(dummy.eigenvec));
% lambda_star = dummy.lambda_star;
% x_star = dummy.x_star;
% sign_FLC = sign(mid(dummy.l1));
%eigenvec = eigenvec / norm(eigenvec);

lambda_star=0.1; init_coord=[-0.060578090252706   0.002256024778798  -0.090978632950962 -0.040896226254254];
T0 = 4.3;

% starting the continuation
[s, x_n] = continuation_Hopf_from_orbit ( lambda_star, fn, f, n_nodes, n_iter, step_size, s, init_coord, T0);
%[s, last_sol] = continuation_Hopf( lambda_star, x_star, f, n_nodes, n_iter, step_size, s, eigenvec, eigenval, sign_FLC);
