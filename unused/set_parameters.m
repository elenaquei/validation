global nu              % weghted ell^1_\nu norm (FLOAT >=1)
global use_intlab      % if intlab will be used or not (BOOL)
global first_run       % distinguish from the first time this function runs (BOOL)
global Display         % if solutions will be plotted or not (BOOL)
%global mu

use_intlab=1;
addpath('../');
addpath('../../');
addpath('./saved elements');
if isempty(first_run) && use_intlab
    addpath('../Intlab_V7.1');
    startintlab;
    first_run=0;
end

% USER'S PARAMETERS  -- you can play
%mu=1.05;       % vanderpol's mu parameter
nu=1.01;        % nu of the l^1_nu space
n_nodes=23;     % number of modes 
delta=3.9e-03;   % length of the first step, from this automatic adaptation
maxiter=60;     % maximum number of iterations of Newton
min_res=2*10^-14; % residual of Newton
set_intlab=0;   % this parameter is used to change from intlab computation to floating point
flag_plot=0;    % plot the numerical solution 
Display =0;     % 0 - doesn't display any value, 1 - display residual of Newton, all the bounds, Imin and Imax
ITER=40;         % maximum number of iterations in the continuation
Max_Comp= 10;   % maximum number of attempted computations of the validation
max_time=Inf; % maximum number of seconds the code will run
delta_decrease = 1.3; % decrease in delta when validation does not work
delta_increase = 1.3; % increase in delta if validation works 5 times in a row


iter_stop=inf;

% Name_system= 'vanderpol_cont';
% Name_system= 'lorenz_cont';
% Name_system = 'lorenz_christian';
% Name_system = 'lorenz_christian_longer';
% Name_system = 'mixingVDP_cont';
Name_system = 'stupidmixing3VDP_cont';
% Name_system = 'fifth_order';
% Name_system = 'lorenz_long';
