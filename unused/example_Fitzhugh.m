% CANNOT FIND A MEANINGFUL ORBIT
%
% let's go towards a Hopf bifurcation in a polynomial ODE
% \dot x = gamma ( x + y - x^3 / 3 + xi )
% \dot y = - ( x - alpha + beta y ) / gamma
% that has a Hopf bifurcation at 
% gamma =1
% beta = 3
% alpha = 0
% xi = 4/3
%     and
% x = 2, y =-2/3
%
% we will continue in xi getting as close as possible to the Hopf
% bifurcation

% definition of some global variables - very important bit
global nu
global RAD_MAX
global Display
global use_intlab
use_intlab = 0;
Display = 0;
RAD_MAX = 10^-5;


try 
    intval(1);
catch
    addpath(genpath('../'))
    startintlab
end

nu = 1.05; 
pert = 0.1;
xi = 4/3-pert;
beta = 3;
n_iter = 10;
s='Fitzhugh';

% DATA OF THE PROBLEM
% number of nodes requested in Fourier
number_of_nodes = 10;

% rhs of the vector field and initial conditions
% initial_coord = [2;-2/3]; unperturbed
initial_coord(2) = -6565;
initial_coord(1) = -3*initial_coord(2);
initial_coord= initial_coord+ 0.1*sqrt(pert);
rhs= @(t,x) -[(x(1)+x(2)-x(1)^3/3+xi);(x(1)+beta* x(2))];
% known approximation of the period of the solution 
approx_period = 2; % even a very rough approximation should do

% string defining the vector field, compatible with the requests of the
% function from_string_to_polynomial_coef 
string_vector_field = ' - dot x1 + l1 x1 + l1 x2 -0.33333 x1^3 +  xi l1 \n   dot x2 + l1  x1 + beta  l1  x2'; 

string_vector_field = strrep(string_vector_field, 'beta', num2str(beta)); % putting in fixed parameters

string_fixed = strrep(string_vector_field, 'xi', num2str(xi)); % for point wise first system, plugging in mu
% continuation vector field, where mu is the second scalar variable
string_cont = strrep(string_vector_field, 'xi', 'l2');


% SETTING UP THE DATA IN RIGHT FORMAT
% computation of the numerical solution with forward integration and built
% in ode solver
%[tout, yout] = ode45(rhs,[0:0.001:approx_period], initial_coord);

% transformation to Xi_vector from time series
%xXi = time_series2Xi_vec(tout,yout,number_of_nodes);

sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;

vector = sqrt(pert)*[sin_four;cos_four];
vector(:,n_nodes+1) = initial_coord;
xXi = Xi_vector(1,vector);

% construction of default scalar equation based on the Xi_vector just
% generated
scalar_equation = fancy_scalar_condition(xXi);

% construction of the vector field from the string
vector_field = from_string_to_polynomial_coef(string_fixed);

% construnction of the full zero-finding problem from the scalar equation
% and the vector field
F = full_problem(scalar_equation, vector_field);

% refinement of the numerical solution 
xXi = Newton_2(xXi,F,30,10^-10);

% defining the continuation problem in pho, starting at pho_null
sol = Xi_vector([xXi.scalar,xi], xXi.vector);
scal_eq = default_scalar_eq(sol);
polynomial = from_string_to_polynomial_coef(string_cont); 
F_not_square = full_problem(scal_eq, polynomial);


% launch the validation
[s, x_n] = continuation ( sol, F_not_square, n_iter, 10^-6, [],s, 10^-10);