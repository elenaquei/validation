% here a certain number of examples are called

% van der pol


%addpath('./saved elements');
%addpath('../Intlab_V7.1');
%startintlab;


%% van der pol from the beginning
init_coord=[-1.288,0.9345];
approx_period=6.3;
initial_mu=1.2;
system='vanderpol';
ITER=36;
n_nodes=30;
nu_input=1.1;
delta_step=10^-3.6;
bool_intlab=1;
bool_display=0;
max_tent=10;
delta_increase=1.2;
delta_decrease=1.3;
bool_invert_dir=1;


[s] = new_thingy (init_coord,initial_mu,approx_period,system,ITER,n_nodes,nu_input,...
    delta_step,bool_intlab,bool_display,max_tent,delta_increase,delta_decrease,...
    bool_invert_dir,[],10^-10);

%% van der pol from mid way
s='/Users/queirolo/Dropbox/amsterdam math/matlab/code/continuation (original)/simulations/vanderpol_cont_validation_86_20.mat';
ITER=5;
[s]=new_thingy(s,ITER);



%% huge lorenz from the beginning, continuation in pho
ITER=1;
nu_input=1.01;
max_Newton_iter=30;
n_nodes=10;
DIM=2; chosed_system=2;
epsilon=10^-8;
[system]=alpha_huge_lorenz(DIM, epsilon);
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

initial_pho=28;
approx_period=possible_T(chosed_system);
init_coord=repmat([possible_x(chosed_system),possible_y(chosed_system),z],1,DIM);

[s] = new_thingy (init_coord,initial_pho,approx_period,system,ITER,n_nodes,nu_input,...
    delta_step,[],[],[],[],[],...
    [],[],10^-12);
load handel
sound(y,Fs)

[s] = new_thingy (init_coord,initial_pho,approx_period,system,ITER,n_nodes,nu_input,...
    delta_step,bool_intlab,bool_display,max_tent,delta_increase,delta_decrease,...
    bool_invert_dir,max_Newton_iter,10^-12);


%% huge lorenz from the beginning, continuation in epsilon
ITER=2;
nu_input=1.01;
max_Newton_iter=30;
n_nodes=30;
DIM=18; 
chosed_system=1;
epsilon=10^-6;
initial_pho=28;
delta_step=10^-3.6;
bool_intlab=1;
bool_display=0;
max_tent=10;
delta_increase=1.2;
delta_decrease=1.3;
bool_invert_dir=1;
% here multiple possibilities

[system]=alpha_epsilon_lorenz(DIM,initial_pho, epsilon);

%[system]=alpha_another_lorenz(DIM,initial_pho);

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

beta=8/3;

[s] = new_thingy (init_coord,beta,approx_period,system,ITER,n_nodes,nu_input,...
    delta_step,bool_intlab,bool_display,max_tent,delta_increase,delta_decrease,...
    bool_invert_dir,max_Newton_iter,10^-14);


%% very long lorenz 
load('saved elements/Periodic_orbit25.mat')
system='lorenz';
ITER=20;
double_c=[flip(conj(c(:,2:end)),2),c];
pho=28;
x=Xi_vector([L,pho],double_c,2,3);
n_nodes=length(c)-1;
nu_input=1.001;
[s] = new_thingy (double_c,pho,L,system,ITER,n_nodes,nu_input,[],[],[],[],[],[],[],[],10^-12);
