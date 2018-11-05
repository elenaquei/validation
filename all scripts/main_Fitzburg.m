% main Fitzhugh
function [s, s_usual_cont]=main_Fitzburg()

global rescaling_saddle 
rescaling_saddle = 500;
global nu
global use_intlab 
global talkative 
global RAD_MAX
talkative = 1;
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;

try 
    intval(1);
catch
    addpath(genpath('./'));
    %addpath(genpath('../../'))
    startintlab;
end

% system data
beta = 3/4;
xi = 5/24;
x = 1/2;
y = -2/3;
eigenval = sqrt(7)/4;
eigenvec = [1/4 *(-3-1i*sqrt(7)),1];
X = [xi,eigenval , x,y,real(eigenvec), imag(eigenvec)]'; % solution
phi = [-4/sqrt(7); -3/sqrt(7)]; % scaling of eigenvalue

s_algebraic = 'Hopf_Fitzhugh'; % where the algebraic solution is saved
algebraicHopf(@system,@derivatives,@second_der,@third_der,2,X,s_algebraic,phi);

% continuation sytem
string_Fitz = '- dot x1 +l1 x1+ l1 x2- 0.3333333 l1 x1^3 + l1 xi \n - dot x2 - l1 x1- beta l1 x2'; % amico fritz
string_Fitz_vars = strrep(string_Fitz, 'beta' , num2str(beta));
string_Fitz_cont = strrep(string_Fitz_vars, 'xi', 'l2'); % setting xi as the second scalar variable


% some elements useful for the computation and the validation
n_nodes = 5; % number of Fourier nodes used: small, since near the Hopf bifurcation is a circle
n_iter = 500; % number of iterations
step_size = 10^-3; % initial step size (then adapted along the validation
s = 'Hopf_Fitz_cont'; % where the solutions are stored

% for Hopf cancel the dependency on the period
vectorfield = strrep(string_Fitz_cont, 'l1' , '');
vectorfield = strrep(vectorfield, 'l2' , 'l1');
% string defining the vector field of the Hopf normal form 

f_Fitz = from_string_to_polynomial_coef(vectorfield); % trasnformation into a vectorfield that can be used

% retrieval of the solution
dummy = load(s_algebraic);

% validated values converted to doubles 
eigenval = conj(mid(dummy.eigenval));
eigenvec = conj(mid(dummy.eigenvec));
lambda_star = dummy.lambda_star;
x_star = dummy.x_star;
sign_FLC = sign(mid(dummy.l1));

% starting the continuation
% [s, last_sol] = continuation_Hopf( lambda_star, x_star, f_Fitz, n_nodes, n_iter, step_size, s, eigenvec, eigenval, sign_FLC);
% equivalent to:
load (s);
last_sol = x_n;

s_usual_cont = 'normal_Fitz_cont';
rescaled_vector = last_sol.scalar(3)*last_sol.vector;
rescaled_vector(:,n_nodes+1) = rescaled_vector(:,n_nodes+1)+last_sol.scalar(4:end)';
back_to_standard_solution = Xi_vector(last_sol.scalar(1:2), rescaled_vector);
f_Fitz_cont = from_string_to_polynomial_coef(string_Fitz_cont);
scal_eq = fancy_scalar_condition(back_to_standard_solution);
F_Fitz_cont = full_problem(scal_eq, f_Fitz_cont);

n_iter_cont=500;
[s_usual_cont, x_n_first] = continuation ( back_to_standard_solution, F_Fitz_cont, 1, -10^-6, [],s_usual_cont, 10^-10);
%xXi = Newton_2(back_to_standard_solution,F_Fitz_cont,30,10^-7);
if x_n_first.scalar(1)<0
    x_n_first.scalar(1) = -x_n_first.scalar(1);
    back_to_standard_solution.scalar(1) = -back_to_standard_solution.scalar(1);
end
figure;plot2(x_n_first,'r'); hold on;plot2(back_to_standard_solution)
load(s_usual_cont)
[s_usual_cont, x_n] = continuation ( x_n, F_Fitz_cont, n_iter_cont, 10^-6, -x_dot_n,s_usual_cont, 10^-10);
if x_n.scalar(1)<0
    x_n.scalar(1) = -x_n.scalar(1);
    back_to_standard_solution.scalar(1) = -back_to_standard_solution.scalar(1);
end
plot2(x_n,'g'); hold on;
end

function f = system(x,xi)
beta = 3/4;
y = x(2);
x = x(1);
f=[(x+y-x^3/3 + xi); 
    -x-beta*y];
end

function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
    DalphaalphaxFn,DxxxFnV] = derivatives(x,xi,v1)
beta = 3/4;
y = x(2);
x = x(1);
DxFn = [1-x^2   1
    -1    -beta];

Dalphafn = [1;0];
if nargin<3
    return
end
DxxFnV = [-2*x*v1(1), 0;0,0];
DalphaxFn = zeros(2);
DalphaxxFn=zeros(2);
DalphaalphaFn=zeros(2,1);
DalphaalphaxFn=zeros(2);
DxxxFnV = zeros(2);
end

function ddf = second_der(x,v1,v2)
beta = 3/4;
gamma=1;
y = x(2);
x = x(1);
ddf = [-2*x*v1(1), 0;0,0]*v2;
end

function dddf = third_der(x,v1,v2,v3)
beta = 3/4;
gamma=1;
dddf = [-2*v1(1)*v2(1)*v3(1);0];
end