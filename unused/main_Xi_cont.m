% MAIN_XI_CONT with the use of the classes Xi_vector, Xi_matrix, coefs
% and coeffs

% STRUCTURE:
% the solutions of the problem are first given by the built-in Matlab ODE
% solver, then it's transformed in Fourier space and the precision is
% augmented by running the Newton method. Finally, the coefficients of the
% radii polynomials and their zeros are computed, and the output is an
% interval [Imin,Imax].

% REMARK : requires Intlab anyway to run

% clear all
global c
c=[];
%clc

%  GLOBAL PARAMETERS  DO NOT TOUCH
global mu
global nu
global use_intlab
global first_run
global Display
use_intlab=1;
addpath('../');
if isempty(first_run) %&& intlab
    addpath('../Intlab_V7.1');
    startintlab;
    first_run=0;
end


% USER'S PARAMETERS  -- you can play
mu=1.2;           % vanderpol's mu parameter ----- continuation parameter
nu=1.06;        % nu of the l^1_nu space
m=60;           % number of modes 
n_nodes=m;
step=10;        % number of continuation steps
delta=0.001;    % length of the first step (FUTURE: automatical adaptation)
maxiter=20;     % maximum number of iterations of Newton
min_res=10^-12; % residual of Newton
set_intlab=0;   % this parameter is used to change from intlab computation to floating point
flag_plot=0;    % plot the numerical solution 
Display =0;     % 0 - doesn't display any value, 1 - display residual of Newton, all th ebounds, Imin and Imax
% Remark: radii polynomials printed anyway

Name_system= 'vanderpol_cont';

fprintf('main_Xi_cont\n\n')
fprintf('\n\nSolving system %s \n\n\n',Name_system);

temp_intlab=use_intlab;
use_intlab=0;

mu0=mu;
mu1=mu+delta;

switch Name_system
    case 'vanderpol_cont'
        
        init_coord=[-1.288,0.9345];
        approx_period=6.3;
        
        % FOURIER COEFS
        disp('van der pol coefficients');
        Name_system='vanderpol';
        interactive_constructor
        
        alpha0=alpha;
        alpha1=alpha;
        coefs_linear0=coefs_linear;
        
        disp('')
        disp('continuation coefficients')
        Name_system='vanderpol_cont';
        interactive_constructor
        
        coefs_linear=coefs_linear0;
        
        alpha0.value{2}=[-2*pi*1i,mu0,-mu0,-1]; 
        % this is due to the dependence
        % on mu, that needs to be computed on the run
        xBarGen0=solve_system(alpha0,m,init_coord,approx_period);
        % numerically solves the problem
        
        
        alpha1.value{2}=[-2*pi*1i,mu1,-mu1,-1]; 
        xBarGen1=solve_system(alpha1,m,init_coord,approx_period);
        % same for x1
        
    otherwise 
        
        fprintf('The system %s was not found, program interrupted.\n',Name_system);
end


if flag_plot
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
    if xBarGen0.size_vector==3
        plot3(xBarGen1,'v');
        hold on
        plot3(xBarGen0,'r');
    else
        plot(xBarGen1,'v');
        hold on
        plot(xBarGen0,'r');
    end
    hold on
end

% NEWTON
[xBar0,yBar0,res0,DH0]=Newton_Xi(xBarGen0,alpha0,coefs_linear,maxiter,min_res);
[xBar1,yBar1,res1,DH1]=Newton_Xi(xBarGen1,alpha1,coefs_linear,maxiter,min_res);


% xBar0=reshape_Xi(xBar0,m*2);
% xBar1=reshape_Xi(xBar1,m*2);
% coefs_linear{2}=cat(3,zeros(1,2,m),cat(3,coefs_linear{2},zeros(1,2,m)));
% m=m*2;

xBar0.size_scalar=2;
xBar0.scalar=[xBar0.scalar,mu0];

xBar1.size_scalar=2;
xBar1.scalar=[xBar1.scalar,mu1];



full_DH0=Function_derivative(xBar0,alpha,coefs_linear,0);
full_DH1=Function_derivative(xBar1,alpha,coefs_linear,0);
% flag 0 to have non-square output

x0_dot=null(full_DH0);
x1_dot=null(full_DH1);

theta_best=0;
for theta=0:0.0001:2*pi
    if norm(x1_dot-exp(1i*theta)*x0_dot)<norm(x1_dot-exp(1i*theta_best)*x0_dot)
        theta_best=theta;
    end
end
x0_dot=exp(1i*theta_best)*x0_dot;


% % other way of computing x0_dot
% mat=[full_DH0;x1_dot.'];
% rhs=[zeros(length(x1_dot)-1,1);1];
% 
% x0_dot2=mat\rhs;


xD_dot=x1_dot-x0_dot;
% E_s (x) = (xBar_s -x) \cdot xs_dot

coefs_linear0=coefs_linear;

coefs_linear0{1}=[x0_dot(1:2).';[0,0]];
coefs_linear0{2}(2,:,:)=reshape(x0_dot(3:end),2,2*m+1);

coefs_linear1=coefs_linear;

coefs_linear1{1}=[x1_dot(1:2).';[0,0]];
coefs_linear1{2}(2,:,:)=reshape(x1_dot(3:end),2,2*m+1);

DH0=[x0_dot.';full_DH0];
DH1=[x1_dot.';full_DH1];


use_intlab=set_intlab;
if use_intlab
    nu=intval(nu);
end


use_intlab=temp_intlab;

%radii_polynomials_cont;
[flag,Imin,Imax] = radii_polynomials_cont(xBar0,xBar1,DH0,DH1,alpha,coefs_linear0, coefs_linear1);

fprintf('\n    The interval found is [%e, %e].\n',Imin,Imax);
