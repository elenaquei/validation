% example_Hopf % working as of 12th September 2017
%
% a code to explain how to validate a Hopf bifurcation 
%
% simple example: the standard Hopf case
%
%   dot x1 = -x2 + x1( lambda - x1^2 -x2^2)
%   dot x2 =  x1 + x2( lambda - x1^2 -x2^2)
% define the system, request the Hopf setting and continuation away from
% the Hopf bifurcation

global nu
global use_intlab
global rescaling_saddle
global refinement_saddle
global talkative
global RAD_MAX
global Display
% setting of some global varibales
Display = 0; % plot of the radii polynomials along the way
talkative = 0; % how much does it talk
use_intlab = 0; % DUMMY, should be changed! 
RAD_MAX = 10^-2; % maximum radius allowed, used for the computation of the Z2 bound
nu = 1.1; % nu of the nu-norm used

rescaling_saddle=1000000;
refinement_saddle = 300;

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
n_nodes = 5; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 10; % number of iterations
step_size = 10^-4; % initial step size (then adapted along the validation
s = 'Hopf_normal'; % where the solutions are stored

vectorfield = '-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3';
%vectorfield = '-dot x1 - x2 + l1 x1 - x1 ^ 3 - x1  x2 ^ 2+x2^4+10*x1^4x2\n- dot x2 + x1 + l1 x2 - x1 ^ 2 x2 - x2 ^ 3+ x1^4';% perturbed for testing purposes
%vectorfield = '-dot x1 + 11x1+37+11x2-3l1-x1^3+9x1^2-31x1-x1x2^2+3x2^2-4x1x2\n-dotx2+13x1+l1x2+2l1-x1^2x2-21x2-29-x2^3-2x1^2_6x1x2-6x2^2'; %translated to (3,-2)

% string defining the vector field of the Hopf normal form 

f = from_string_to_polynomial_coef(vectorfield); % trasnformation into a vectorfield that can be used

% definition of the solution
x_star =[0,0];
%x_star =[3,-2];
lambda_star = 0;

% starting the continuation
eigenvec = [1/2*1i,1/2]+ rand*0.001;
eigenval = 1i;

[s, last_sol] = continuation_Hopf( lambda_star, x_star, f, n_nodes, n_iter, step_size, s, eigenvec, eigenval);
