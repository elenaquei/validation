% high dimension van der Pol with new mixing (i.e. non-diagonal problem)
% x dot = y_1
% y_i dot = mu (1-x^2)*(y_i+y_(i+1))/2-x
global nu
global RAD_MAX
global Display
global use_intlab
use_intlab = 0;
Display = 1;
RAD_MAX = 10^-8;

try 
    intval(1);
catch
    addpath(genpath('../'))
    startintlab
end

nu = 1.05; 

% DATA OF THE PROBLEM
% number of nodes requested in Fourier
number_of_nodes = 60;
DIM= 3;

% rhs of the vector field and initial conditions
initial_coord = [2,0];
rhs= @(t,x) [x(2);(1-x(1).^2)*x(2)-x(1)];
% known approximation of the period of the solution 
approx_period = 6.8; % even a very rough approximation should do

% string defining the vector field, compatible with the requests of the
% function from_string_to_polynomial_coef 
string_vector_field = 'dot x1 - l1 x2 \n dot x2 - l1 x2 + l1 x1^2 x2 + l1 x1'; 


% SETTING UP THE DATA IN RIGHT FORMAT
% computation of the numerical solution with forward integration and built
% in ode solver
[tout, yout] = ode45(rhs,[0:0.001:approx_period], initial_coord);

% transformation to Xi_vector from time series
xXi = time_series2Xi_vec(tout,yout,number_of_nodes);

% construction of the vector field from the string
vector_field = coupling_VDP(DIM);

% make it big!
xXi_big=xXi;
xXi_big.size_vector = 1+DIM;
vector2 = repmat(xXi.vector(2,:),DIM,1);
xXi_big.vector = [xXi.vector(1,:);vector2];

% construction of default scalar equation based on the Xi_vector just
% generated
scalar_equation = default_scalar_eq(xXi_big,2,0.4);
scalar_equation = fancy_scalar_condition(xXi_big, scalar_equation, 1);

% construnction of the full zero-finding problem from the scalar equation
% and the vector field
F = full_problem(scalar_equation, vector_field);

% refinement of the numerical solution 
xXi = Newton_2(xXi_big,F,30,10^-10);

% VALIDATION
success = validation_orbit(F, xXi);
