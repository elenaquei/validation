% 31 August 2017
%
% non-coupled works fine
% coupled works fine BUT up to DIM<5 (at 5 the continuation fails, with
% this parameters) 2 and 3 validated (far from 19...)

function s = coupled_lorenz(DIM)
% function s = coupled_lorenz(DIM)
% 
% validated cntinuation of DIM coupled Lorenz systems, starting from a
% well-known orbit
    
% standard definition of necessary elements
global azabaza
global nu
global use_intlab
global talkative
global RAD_MAX
global Display
Display = 1;
talkative = 5;
use_intlab = 0;
RAD_MAX = 10^-2;

if isempty(azabaza)
    addpath(genpath('../'))
    startintlab
    azabaza =1;
end

% problem dependent
nu = 1.1;
pho_null = 28;
n_nodes = 50;
step_size = 10^-4;
sigma = 10;
beta = 8/3;
pho = 28;

chosed_system=2;


% construct the numerical solution with forward integration from known
% initial conditions
possible_x=[-1.3763610682134e+01; 
    -1.2595115397689e+01; 
    -1.4426408025035e+01;
    -1.3056930146345e+01];
possible_y=[-1.9578751942452e+01
    -1.6970525307084e+01
    -2.1111230056994e+01
    -1.7987214049281e+01];
possible_T=[1.5586522107162e+00
    2.3059072639398e+00
    2.3059072639398e+00
    3.8202541634368e+00];
z=27;

approx_period=possible_T(chosed_system);
init_coord=repmat([possible_x(chosed_system),possible_y(chosed_system),z],1,1);

one_dim_coord = [possible_x(chosed_system),possible_y(chosed_system),z];
f=@(t,x)[sigma*(x(2) - x(1)); x(1)*(pho-x(3)) - x(2); x(1)*x(2)- beta*x(3)];

% forward integration
[tout, yout] = ode45(f,[0,approx_period],one_dim_coord);

% transformation to Xi_vector from time series
xXi = time_series2Xi_vec(tout,yout,n_nodes);

%% standard Lorenz continuation

% definition of the vector field in the form of a string
string_lorenz = '- dot x1 + sigma l1 x2 - sigma l1 x1 \n - dot x2 + pho l1 x1 - l1 x1 x3 - l1 x2 \n - dot x3 + l1 x1 x2 - beta l1 x3'; % general lorenz
string_lorenz_vars = strrep(string_lorenz, 'sigma' , num2str(sigma)); % plugging in sigma
string_lorenz_vars = strrep(string_lorenz_vars, 'beta' , num2str(beta)); % plugging in beta
string_lorenz_pho = strrep(string_lorenz_vars, 'pho', num2str(pho_null)); % for point wise first system, plugging in pho
string_lorenz_cont = strrep(string_lorenz_vars, 'pho', 'l2'); % setting pho as the second scalar variable

% constructing solution to the problem
sol = xXi;%
scal_eq = default_scalar_eq(sol);

% constructing the ODEs systems - the point and the continuous
% fixed pho
polynomial_fix = from_string_to_polynomial_coef(string_lorenz_pho);
F_fix = full_problem(scal_eq, polynomial_fix);

% NEWTON
[sol2,yBar,res,DFm,RES] =Newton_2(sol,F_fix,30,10^-7);

% continuation

% solution at pho_null
sol = Xi_vector([sol2.scalar,pho_null], sol2.vector);

% non-square problem
scal_eq = default_scalar_eq(sol);
polynomial = from_string_to_polynomial_coef(string_lorenz_cont); 
F_not_square = full_problem(scal_eq, polynomial);

string_coupled = create_coupled_system(DIM);
sol.size_vector = sol.size_vector*DIM;
sol.vector = repmat(sol.vector,DIM,1);
scal_eq = default_scalar_eq(sol);
polynomial = from_string_to_polynomial_coef(string_coupled); 
F_coupled = full_problem(scal_eq, polynomial);


% test on a point, the continuation does not work




% call the continuation function
n_iter =1;
s = sprintf('saved elements/coupled_lorenz_DIM%i', DIM);
[s, x_n] = continuation ( sol, F_coupled, n_iter, step_size, [],s, 10^-10);
end


function s = create_coupled_system(dim)
% function s = create_coupled_system(dim)
% 
% INPUT
% dim    integer, number of coupled systems
% OUTPUT
% s      string, dim Lorenz coupled systems
%
% continuation in pho

sigma = 10;
beta = 8/3;

base =  '- dot xa + sigma l1 xb - sigma l1 xe \n - dot xb + pho l1 xa - l1 xa xc - l1 xb \n - dot xc + l1 xa xb - beta l1 xc\n'; 
base = strrep(base, 'sigma' , num2str(sigma)); % plugging in sigma
base = strrep(base, 'beta' , num2str(beta)); % plugging in beta
base = strrep(base, 'pho', 'l2'); 

s ='';
i=1;
step = strrep(base,'a',num2str((i-1)*3+1));
step = strrep(step,'b',num2str((i-1)*3+2));
step = strrep(step,'c',num2str((i-1)*3+3));
step = strrep(step,'e',num2str((dim-1)*3+1));
s = strcat(s,step);

for i =2:dim
    step = strrep(base,'a',num2str((i-1)*3+1));
    step = strrep(step,'b',num2str((i-1)*3+2));
    step = strrep(step,'c',num2str((i-1)*3+3));
    step = strrep(step,'e',num2str((i-1)*3-2));
    s = strcat(s,step);
end




end
