% some examples to run the code and check it out
% 31 August 2017: runs

global nu
global use_intlab
global talkative
global RAD_MAX
global Display
Display = 0;
talkative = 1;
use_intlab = 0;
RAD_MAX = 10^-2;

try 
    intval(1);
catch
    addpath(genpath('../'))
    startintlab
end


%% standard Lorenz continuation

nu = 1.1;
pho_null = 28;
n_nodes = 50;
step_size = 10^-2;
sigma = 10;
beta = 8/3; 

string_lorenz = '- dot x1 + sigma l1 x2 - sigma l1 x1 \n - dot x2 + pho l1 x1 - l1 x1 x3 - l1 x2 \n - dot x3 + l1 x1 x2 - beta l1 x3'; % general lorenz
string_lorenz_vars = strrep(string_lorenz, 'sigma' , num2str(sigma)); % plugging in sigma
string_lorenz_vars = strrep(string_lorenz_vars, 'beta' , num2str(beta)); % plugging in beta
string_lorenz_pho = strrep(string_lorenz_vars, 'pho', num2str(pho_null)); % for point wise first system, plugging in pho
string_lorenz_cont = strrep(string_lorenz_vars, 'pho', 'l2'); % setting pho as the second scalar variable

% loading approximate solution 

load('huge_lorenz1_validation_1_40','x0');
x0.scalar(1) = x0.scalar(1)/(2*pi);
x0 = reshape_Xi(x0,n_nodes);

% constructing solution to the problem
sol = x0;% Xi_vector(1/L, [x1(2:end);x2(2:end); x3(2:end)]);
sol.scalar = sol.scalar(1);
sol.size_scalar = 1;
scal_eq = default_scalar_eq(sol);

% constructing the ODEs systems - the point and the continuous
% fixed pho
polynomial_fix = from_string_to_polynomial_coef(string_lorenz_pho);
F_fix = full_problem(scal_eq, polynomial_fix);

% NEWTON
[sol2,yBar,res,DFm,RES] =Newton_2(sol,F_fix,30,10^-7);

% validation poitwise
DF =  derivative(F_fix,sol2,0);
DF_mat = derivative_to_matrix(DF);
A  = inv(DF_mat);


use_intlab = 1;

Y_vector = Y_bound_new(A,sol2,F_fix);
Z0_vector=Z0_bound(DF_mat,A,sol2);
Z1_vector=Z1_bound_new(A,sol2,F_fix);
Z2_vector= Z2_bound_new(A,sol2,F_fix);
[Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);

% check
if Imax>RAD_MAX
    Imax = RAD_MAX;
end

if talkative
    fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
end

% use_intlab = 0;


% continuation
% in pho, starting at pho_null
sol = Xi_vector([sol2.scalar,pho_null], sol2.vector);

scal_eq = fancy_scalar_condition(sol);
polynomial = from_string_to_polynomial_coef(string_lorenz_cont); 

F_not_square = full_problem(scal_eq, polynomial);

n_iter =10000;
s = './saved elements/lorenz_continuation';
[s, x_n] = continuation ( sol, F_not_square, n_iter, step_size,[],s,10^-10, 0,1, 0);

% still not thaaaat cool
% (no fold in pho yet)
