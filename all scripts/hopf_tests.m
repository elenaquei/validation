function hopf_tests
% example_Hopf % working as of 12th September 2017
%
% a code to explain how to validate a Hopf bifurcation 
%
% the standard Hopf case
%
%   dot x1 = -x2 + x1( lambda - x1^2 -x2^2)
%   dot x2 =  x1 + x2( lambda - x1^2 -x2^2)
% define the system, request the Hopf setting and continuation away from
% the Hopf bifurcation

global nu
global use_intlab
global rescaling_saddle
global refinement_saddle
%global talkative
global RAD_MAX
global Display
% setting of some global varibales
Display = 0; % plot of the radii polynomials along the way
%talkative = 10; % how much does it talk
use_intlab = 0; % DUMMY, should be changed! 
RAD_MAX = 10^-2; % maximum radius allowed, used for the computation of the Z2 bound
nu = 1.1; % nu of the nu-norm used

rescaling_saddle=10^4;
refinement_saddle = 700;

% if first run, add Intlab and start it (take care, starting Intlab cancels
% all break points)
try 
    intval(1);
catch
    addpath(genpath('./'))
    addpath(genpath('../'))
    startintlab
end

% some elements useful for the computation and the validation
n_nodes = 7; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 10; % number of iterations
step_size = 10^-4; % initial step size (then adapted along the validation
s = 'Hopf_small_test'; % where the solutions are stored

vectorfield = '-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3';

% string defining the vector field of the Hopf normal form 

f = from_string_to_polynomial_coef(vectorfield); % trasnformation into a vectorfield that can be used

% definition of the solution
x_star =[0,0];
%x_star =[3,-2];
lambda_star = 0;

% starting the continuation
eigenvec = [1/2*1i,1/2]+ rand*0.001;
eigenval = 1i;

big_Hopf = construct_big_hopf( ...
       lambda_star, x_star, f, n_nodes, eigenvec, eigenval);

[~, x_n] = continuation_Hopf( ...
    lambda_star, x_star, f, n_nodes, 1, step_size, s, eigenvec, eigenval);
save(s,'x_n','big_Hopf','-append')
load(s)
%load Hopf_normal


s_other = 'hopf_continuing';
step_size = 10^-4;
save_x_n = 'hopf_solutions';
num_simul = 100;

% plot_transform(x_n)
% hold on

run_simulation = 1;
if run_simulation
    x_n_cell = cell(num_simul,1);
    for i = 1:num_simul
        [s_other, x_n] = continuation ( x_n, big_Hopf, n_iter, step_size, x_dot_n, s_other, 10^-6);
        % plot_transform(x_n)
        x_n_cell{i} = x_n;
        disp(i)
    end
    save(save_x_n,'x_n_cell');
end

load(save_x_n)

vec_plots = [1:6:num_simul/3,floor(num_simul/3)+2:3:num_simul];%, num_simul/2+4:4:num_simul] ;
vec_plots = vec_plots(1:2:length(vec_plots));
for i = vec_plots
    plot_transform(x_n_cell{i})
    hold on
end

return

end

% order of prarameters
% [ 1/T, lambda, a, x ]

function  plot_transform(x)

y = x;

y.vector = x.scalar(3)*x.vector;

y.vector(:,x.nodes+1)=x.scalar(4:end).';

y_time_series = time_series(y);

lambda = x.scalar(2) + 0*y_time_series(:,1);

plot3(lambda, y_time_series(:,1), y_time_series(:,2),'LineWidth',2);
hold on
plot3(x.scalar(2), x.scalar(4), x.scalar(5),'.','LineWidth',2);
end





