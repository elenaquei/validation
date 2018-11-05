function [x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
    max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
    Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,bool_plot,Max_comp,Name_system,delta_increase,delta_decrease]...
    = load_elements(X0,var0)%,period0,system,n_nodes_new)
% function [x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
%     max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
%     Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,bool_plot,Max_comp,Name_system,delta_increase,delta_decrease]...
%     = load_elements(X0,var0)%,period0,system,n_nodes_new)

global nu              % weghted ell^1_\nu norm (FLOAT >=1)
global use_intlab      % if intlab will be used or not (BOOL)
global first_run       % distinguish from the first time this function runs (BOOL)
global Display         % if solutions will be plotted or not (BOOL)
global talkative       % if talkative==0 no output
                       % if talkative==1 some output
                       % if talkative==2 more output

ITER_new=var0;

load(X0);
ITER=ITER_new+ITER;

delta_step=delta_vec(iter);
delta_vec=[delta_vec(1:iter);zeros(ITER_new,1)];
Imin_vec=[Imin_vec(1:iter);zeros(ITER_new,1)];
Imax_vec=[Imax_vec(1:iter);zeros(ITER_new,1)];
xDelta=[xDelta(1:iter);zeros(ITER_new,1)];
num_nodes_vec=[num_nodes_vec(1:iter);zeros(ITER_new,1)];
tot_size=x0.size_scalar+x0.size_vector;
Y0_vec=[Y0_vec(1:iter,:);zeros(ITER_new,tot_size)];
Z0_vec=[Z0_vec(1:iter,:);zeros(ITER_new,tot_size)];
Z1_vec=[Z1_vec(1:iter,:);zeros(ITER_new,tot_size)];
Z2_vec=[Z2_vec(1:iter,:);zeros(ITER_new,tot_size)];
normX1_vec=[normX1_vec(1:iter+1,:);zeros(ITER_new,tot_size)];
pho_vec=[pho_vec(1:iter+1);zeros(ITER_new,1)];
first_iter=iter;
%clear ITER_new
if ~exist('max_Newton_iter','var')
    max_Newton_iter=60;
end
if ~exist('res_Newton','var')
    res_Newton=10^-14;
end

if ~exist('bool_plot','var')|| isempty(bool_plot)
    bool_plot=0;
end
if ~exist('end_time','var')
    end_time=0;
end
if ~exist('Max_comp','var')
    Max_comp=10;
end
if ~exist('DH0','var')
    DH0=Function_derivative(x0,alpha_coef,coefs_linear0);
end
if isempty(talkative)
    if length(X0)*(2*n_nodes+1)>200
        talkative=2;
    else
        talkative=1;
    end
end
if talkative>0
    fprintf('Starting the code, time %s', datetime(now,13))
end
return
end