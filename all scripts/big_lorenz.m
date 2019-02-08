function [s, xn] = big_lorenz(DIM, num_iter)
% lorenz problem
% b = 8/3, sigma = 10, p = 28
%
% x_i dot = sigma( y_i - x_i)
% y_i dot = x_i ( rho + eps*i - z_i ) - y_(i+1)
% z_i dot = x_i y _i - beta z_i

global first_run
global nu
global use_intlab 
global talkative 
global RAD_MAX
talkative = 1;
use_intlab = 0;
nu = 1.01;
RAD_MAX = 10^-4;

if isempty(first_run)
    addpath(genpath('../'))
    startintlab
    first_run =1;
end

n_nodes= 65;
sigma=10;beta=8/3;
pho_null =28;

string_lorenz = '- dot x1 + sigma l1 x2 - sigma l1 x1 \n - dot x2 + pho l1 x1 - l1 x1 x3 - l1 x2 \n - dot x3 + l1 x1 x2 - beta l1 x3'; % general lorenz
string_lorenz_vars = strrep(string_lorenz, 'sigma' , num2str(sigma)); % plugging in sigma
string_lorenz_vars = strrep(string_lorenz_vars, 'beta' , num2str(beta)); % plugging in beta
string_lorenz_cont = strrep(string_lorenz_vars, 'pho', 'l2'); % setting pho as the second scalar variable


chosed_system=2;

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
sigma = 10;
beta = 8/3;
pho = 28;

one_dim_coord = [possible_x(chosed_system),possible_y(chosed_system),z];
f=@(t,x)[sigma*(x(2) - x(1)); x(1)*(pho-x(3)) - x(2); x(1)*x(2)- beta*x(3)];

[tout, yout] = ode45(f,[0,approx_period],one_dim_coord);
%plot3(yout(:,1),yout(:,2),yout(:,3))

size_vec = 3;
size_scalar =2;

x=0*yout;
for i=1:size_vec
    x(:,i)=1/(size(yout,1))*fft(yout(:,i));
end
xBar=cell(size_vec+size_scalar,1);
xBar{1}=approx_period/(2*pi);%time;
xBar{2}=0;
m=n_nodes;

for i=1:size_vec
    xBar{i+size_scalar}=fftshift([x(1:m+1,i);x(end-m+1:end,i)]);
    xBar{i+size_scalar}(1:m)=conj(flip(xBar{i+size_scalar}(m+2:end)));
    xBar{i+size_scalar}(m+1)=real(xBar{i+size_scalar}(m+1));
end

scal=[approx_period^-1,pho];
vec=[xBar{1+size_scalar}];
for i=2:size_vec
    vec=[vec,xBar{i+size_scalar}];
end
xXi=Xi_vector(scal,vec,size_scalar,size_vec,n_nodes);


scal_eq = default_scalar_eq(xXi);
scal_eq = fancy_scalar_condition(xXi,scal_eq,1);

% constructing the ODEs systems - the point and the continuous
% fixed pho
polynomial_cont = from_string_to_polynomial_coef(string_lorenz_cont);
F_lor_non_square = full_problem(scal_eq, polynomial_cont);

% norm(apply(F_lor,x0))

% set the two scalar equations

% pho constant - soon removed
F_lor = F_lor_non_square;
DF=derivative_to_matrix(derivative(F_lor,xXi,0));
x_dot0=kernel(DF);
coefs_linear_lor = cell(1,3);
coefs_linear_lor{1} = [0,1];
coefs_linear_lor{2} = 0*shiftdim(reshape(-x_dot0(xXi.size_scalar+1:end),2*n_nodes+1,xXi.size_vector).',-1);
coefs_linear_lor{3} = -pho_null; 
F_lor.scalar_equations = change_lin_coef(...
    F_lor.scalar_equations,coefs_linear_lor,2);

xXi = Newton_2(xXi,F_lor,[],10^-11);

% make it big!
xXi_big=xXi;
xXi_big.scalar(2) = 0;
xXi_big.size_vector = DIM*3;
xXi_big.vector = repmat(xXi.vector,DIM,1);

polynomial_VF = vector_field_lorenz_mixing(DIM,[]);
scal_eq = default_scalar_eq(xXi_big);
F_lor_not_square = full_problem(scal_eq, polynomial_VF);

use_intlab=0;
s =strcat(pwd,strcat('/coupled_lorenz',num2str(DIM)));
[s, xn] = continuation(xXi_big,F_lor_not_square,num_iter, 10^-4, [],s ,10^-11);
% the results will be stored in s

