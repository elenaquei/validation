global first_run
global nu
global use_intlab 
global talkative 
global RAD_MAX
talkative = 5;
use_intlab = 0;
nu = 1.01;
RAD_MAX = 10^-2;

if isempty(first_run)
    addpath(genpath('../'))
    addpath(genpath('../../'))
    startintlab
    first_run =1;
end

% interesting values
d = 1;
f = 1;
zeta = 1;
epsilon = 0.001;

n_nodes = 80;

sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;
y_vec = epsilon * [sin_four; cos_four;sin_four; cos_four];
x0 = Xi_vector(2*pi,y_vec);

string_vf = 'dot x1 - dinv l1  x4 - zeta dinv l1 x2 + dinv l1 x2^3 + dinv l1 x2 x4^2 \n dot x2 - l1 x1 \n \dot x3 + dinv l1 x2 - zeta dinv l1 x4 + dinv l1 x2^2 x4 + dinv l1 x4^3 - dinv f l1 \n dot x4 - l1 x3';

string_vf = strrep(string_vf, 'zeta dinv', num2str(zeta/d));
string_vf = strrep(string_vf, 'dinv f', num2str(f/d));
string_vf = strrep(string_vf, 'dinv', num2str(1/d));

f = from_string_to_polynomial_coef(string_vf);

scalar_eqs_fixed = fancy_scalar_condition(x0);
F_fixed = full_problem(scalar_eqs_fixed,f);
x0_N = Newton_2(x0,F_fixed);  % poor choice, probably.

[flag,Imin,Imax]=radii_polynomials(x0_N,F_fixed);
