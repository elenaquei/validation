% doesn't work and double of examples_lorenz (31 August 2017)


% lorenz problem
% b = 8/3, sigma = 10, p starts at 28
%
% dx = sigma * y - sigma * x)
% dy = p * x - x * z - y
% dz = x * y - b * z


global first_run
global nu
global use_intlab 
global talkative 
global RAD_MAX
talkative = 2;
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;
n_iter = 3;

if isempty(first_run)
    addpath(genpath('../'))
    startintlab
    first_run =1;
end

n_nodes= 60;
DIM = 1;
sigma=10;beta=8/3;

%F_lor = new_alpha_huge_lorenz(DIM,[],n_nodes);
F_lor = new_alpha_lorenz([],n_nodes);


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
init_coord=repmat([possible_x(chosed_system),possible_y(chosed_system),z],1,DIM);
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
xXi=Xi_vector(scal,repmat(vec,1,DIM),size_scalar,DIM*size_vec,n_nodes);

load('huge_lorenz1_validation_1_40','x0');
x0.scalar(1) = x0.scalar(1)/(2*pi);
x0 = reshape_Xi(x0,n_nodes);

% norm(apply(F_lor,x0))

% set the two scalar equations

% pho constant - soon removed
DF=derivative_to_matrix(derivative(F_lor,xXi,0));
x_dot0=kernel(DF);
coefs_linear_lor = cell(1,3);
coefs_linear_lor{1} = [0,1];
coefs_linear_lor{2} = 0*shiftdim(reshape(-x_dot0(xXi.size_scalar+1:end),2*n_nodes+1,xXi.size_vector).',-1);
coefs_linear_lor{3} = -xXi.scalar(2); 
F_lor.scalar_equations = change_lin_coef(...
    F_lor.scalar_equations,coefs_linear_lor,2);

% z(0) constant (w.r.t. xXi) - kept
coefs_linear_lor = cell(1,3);
coefs_linear_lor{1} = 0*[1,1];
coefs_linear_lor{2} = 0*shiftdim(reshape(-x_dot0(xXi.size_scalar+1:end),2*n_nodes+1,xXi.size_vector).',-1);
coefs_linear_lor{2}(1,3,:)=1;
coefs_linear_lor{3} = real(-sum((xXi.vector(3,:))));
F_lor.scalar_equations = change_lin_coef(...
    F_lor.scalar_equations,coefs_linear_lor,1);


xXi = Newton_2(x0,F_lor,[],10^-11);

F_lor.scalar_equations = remove_lin_coef(F_lor.scalar_equations,xXi.size_scalar);
s = '/Users/queirolo/Desktop/simulations/lorenz.mat';
[s, xn] = continuation(xXi,F_lor,n_iter, 10^-3, [],s ,10^-11);
% the results will be stored in s

