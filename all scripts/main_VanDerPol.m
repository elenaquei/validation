% Basic Van Der Pol system, 1 validation
% as of 26 October 2018: it runs

global azabaza
global nu
global use_intlab
global talkative
global RAD_MAX
global Display
Display = 1;
talkative = 1;
use_intlab = 0;
RAD_MAX = 10^-2;

if isempty(azabaza)
    addpath(genpath('../'))
    startintlab
    azabaza =1;
end


%% standerad VDP
% Van der Pol 


mu_vector = 1;
nu_vector = 1.001;
n_nodes_vector = 50;
power_vector = 1;

string_van_der_pol = '- dot x1 + l1 x2 \n - dot x2 + mu l1 x2 - mu l1 x1^2 x2 - l1 x1'; % general van der pol
n_nodes = 50;
% construction of approximate solution (taking the circle every time, easy)
sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
x1 = sin_four;
x2 = cos_four;

% constructing the problem, vector field and simple phase condition
sol = Xi_vector(1, [x1;x2]);
sol2 = sol;

nu = nu_vector(1);
mu = mu_vector(1);


string_van_der_pol_mu = strrep(string_van_der_pol, 'mu', num2str(mu)); % plugging in mu
polynomial = from_string_to_polynomial_coef(string_van_der_pol_mu);

use_intlab = 0;

sol2 = reshape_Xi(sol2,n_nodes);
scal_eq = default_scalar_eq(sol2);
F = full_problem(scal_eq, polynomial);

% NEWTON
try
    [sol2,yBar,res,DFm,RES] =Newton_2(sol2,F,30,10^-7);
    
catch
    fprintf('Newton failed \n')
    break
end

% validation
DF =  derivative(F,sol2,0);
DF_mat = derivative_to_matrix(DF);
A  = inv(DF_mat);

use_intlab = 1;

try
    T1 = cputime;
    Y_vector = Y_bound_new(A,sol2,F);
    Z0_vector=Z0_bound(DF_mat,A,sol2);
    Z1_vector=Z1_bound_new(A,sol2,F);
    Z2_vector= Z2_bound_new(A,sol2,F);
    [Imin,Imax]=find_negative(Z2_vector,Z1_vector,Z0_vector,Y_vector);
    
    % check
    if Imax>RAD_MAX
        Imax = RAD_MAX;
    end
    
    if talkative
        fprintf('\n The interval is [ %e, %e ].\n\n',Imin,Imax);
    end
    
catch
end

return