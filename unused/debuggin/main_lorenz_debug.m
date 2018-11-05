% lorenz problem
% b = 8/3, sigma = 10, p starts at 28
%
% dx = sigma * y - sigma * x)
% dy = p * x - x * z - y
% dz = x * y - b * z


global azabaza
global nu
global use_intlab 
global talkative 
global RAD_MAX
talkative = 1;
use_intlab = 0;
nu = 1.1;
RAD_MAX = 10^-4;

if isempty(azabaza)
    addpath(genpath('../'))
    startintlab
    azabaza =1;
end

n_nodes= 40;
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

norm(apply(F_lor,x0))

%debugging
one_Xi = zeros(1,2*n_nodes+1);
one_Xi(n_nodes+1) =1;
Conv = @(x,y) conv(x,y,'same');
der = 1i*(-x0.nodes:x0.nodes);
compare = @(x,y,z,per,pho) [der.*x-per*sigma*(y-x);der.*y-per*pho*x+per*Conv(x,z)+per*y;der.*z-per*Conv(x,y)+per*beta*z];
compare_Xi = @(Xi) compare(Xi.vector(1,:),Xi.vector(2,:),Xi.vector(3,:),Xi.scalar(1),Xi.scalar(2));
%xXi = x0;
%xXi.vector = compare_Xi(x0);

% set the two scalar equations
DF=derivative_to_matrix(derivative(F_lor,xXi,0));
x_dot0=kernel(DF);
coefs_linear_lor = cell(1,3);
coefs_linear_lor{1} = [0,1];%-x_dot0(1:xXi.size_scalar).';
coefs_linear_lor{2} = 0*shiftdim(reshape(-x_dot0(xXi.size_scalar+1:end),2*n_nodes+1,xXi.size_vector).',-1);
coefs_linear_lor{3} = -xXi.scalar(2); %0*Xi_vec2vec(xXi).'*x_dot0;
F_lor.scalar_equations = change_lin_coef(...
    F_lor.scalar_equations,coefs_linear_lor,2);

coefs_linear_lor = cell(1,3);
coefs_linear_lor{1} = 0*[1,1];%-x_dot0(1:xXi.size_scalar).';
coefs_linear_lor{2} = 0*shiftdim(reshape(-x_dot0(xXi.size_scalar+1:end),2*n_nodes+1,xXi.size_vector).',-1);
coefs_linear_lor{2}(1,3,:)=1;
coefs_linear_lor{3} = -sum((xXi.vector(3,:)));
F_lor.scalar_equations = change_lin_coef(...
    F_lor.scalar_equations,coefs_linear_lor,1);



coefs_linear_lor{3} = -sum((x0.vector(3,:)));
F_lor.scalar_equations = change_lin_coef(...
    F_lor.scalar_equations,coefs_linear_lor,1);
%x0.scalar(2) =0;
%x0.scalar(1) = 2.3/(2*pi);

norm(apply(F_lor,x0))
DF_der = derivative(F_lor,x0,0);
DF_CODE = derivative_to_matrix(DF_der);
DF_fin_dif = 0*DF_CODE;
h = 0.00000001;
e = eye(length(DF_fin_dif));
for i=1:size(DF_fin_dif)
    x_h = vec2Xi_vec(Xi_vec2vec(x0)+h*e(:,i),x0);
    DF_fin_dif(:,i) = Xi_vec2vec(apply(F_lor,x0)-apply(F_lor,x_h))/(-h);
end
pad = zeros(1,n_nodes);
x = xXi.vector(1,:);
y = xXi.vector(2,:);
z = xXi.vector(3,:);
K = -n_nodes:n_nodes;
inverse_conj= @(vec) conj(vec(end:-1:1));
second_half = @(vec) vec(((length(vec)+1)/2):end);
first_half = @(vec) vec(((length(vec)+1)/2:-1:1));
der=@(vec) conj(toeplitz(([first_half(vec),pad]),([second_half(vec),pad])));
pho = xXi.scalar(2);
per = xXi.scalar(1);

DF_ex = [diag(1i*K)+der(per*sigma*one_Xi), der(-per*sigma*one_Xi), der(0*one_Xi);
    der(-per*pho*one_Xi+per*z), diag(1i*K)+der(per*one_Xi), der(per*x);
    der(-per*y), der(-per*x),diag(1i*K)+ der(per*beta*one_Xi)];

xXi = Newton_2(x0,F_lor); % WORKS, SILLY ME!


return
% test problem
value = cell(2,1);
dot = cell(2,1);
power_vec = cell(2,1);
power_scal = cell(2,1);
power_scal{1} = [0,1];
power_scal{2} = [0,1];
value{1} = [-1,1];
value{2} = [-1,-1];
dot{1}= cell(2,1);
power_vec{1}= cell(2,1);
dot{2}= cell(2,1);
power_vec{2}= cell(2,1);
dot{1}{1} = [1,0].';
dot{1}{2} = [0,0].';
power_vec{1}{1} = [0,0].';
power_vec{1}{2} = [0,1].';
dot{2}{1} = [0,1].';
dot{2}{2} = [0,0].';
power_vec{2}{1} = [0,0].';
power_vec{2}{2} = [1,0].';
F_pol = polynomial_coefs(1,2,2,[2,2],value,power_scal,power_vec,dot);
coefs_linear0{2} = 1+coefs_linear0{2}(1,1:2,:);
coefs_linear0{1} = 0;
coefs_linear0{3} = -1;

pol = polynomial_coefs(1,2,0,zeros(0),cell(0),cell(0),cell(0),cell(0));

F_scal = scalar_eq(1,0,1,2,coefs_linear0,pol);

F_test = full_problem(F_scal,F_pol);

size_vec = 2;
size_scalar =1;

x=[sin(0:2*pi/(2*n_nodes+1):2*pi);cos(0:2*pi/(2*n_nodes+1):2*pi)].';
for i=1:size_vec
    x(:,i)=1/(size(x,1))*fft(x(:,i));
end


xBar=cell(size_vec+size_scalar,1);
xBar{1}=1;%time;
xBar{2}=0;
m=n_nodes;

for i=1:size_vec
    xBar{i+size_scalar}=fftshift([x(1:m+1,i);x(end-m+1:end,i)]);
    xBar{i+size_scalar}(1:m)=conj(flip(xBar{i+size_scalar}(m+2:end)));
    xBar{i+size_scalar}(m+1)=real(xBar{i+size_scalar}(m+1));
end

scal=1;
vec=[xBar{1+size_scalar}];
for i=2:size_vec
    vec=[vec,xBar{i+size_scalar}];
end
sin_four = zeros(1,2*n_nodes+1);
cos_four = sin_four;
cos_four(n_nodes) = 1/2; cos_four(n_nodes+2) = 1/2;
sin_four(n_nodes) = 1i/2; sin_four(n_nodes+2) = -1i/2;

vec2 = [sin_four; cos_four];

xXi=Xi_vector(1,vec2,size_scalar,size_vec,n_nodes);


norm(apply(F_test,xXi))
