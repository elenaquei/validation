function [x0,DH0,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
    max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
    Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec]=compute_elements(X0,var0,period0,system,ITER,n_nodes,...
    bool_invert_dir,max_Newton_iter,res_Newton,bool_plot)
% function [x0,DH0,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
%     max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
%     Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec]=compute_elements(X0,var0,period0,system,ITER,n_nodes,...
%     bool_invert_dir,max_Newton_iter,res_Newton,bool_plot)

global nu              % weghted ell^1_\nu norm (FLOAT >=1)
global use_intlab      % if intlab will be used or not (BOOL)
global first_run       % distinguish from the first time this function runs (BOOL)
global Display         % if solutions will be plotted or not (BOOL)


end_time=0;
%delta=delta_step;
%maxiter=max_Newton_iter;
%min_res=res_Newton;
if ischar(system)
    Name_system=system;
    % interactive_constructor is a script that recalls the already
    % stored coefficients of the vector field or interactively
    % construct a new vector field
    interactive_constructor
    
    alpha0=alpha_coef;
    
    coefs_linear0=coefs_linear;
    
    %xBarGen0=solve_system(alpha0,n_nodes,init_coord,approx_period);
    
    disp('')
    disp('continuation coefficients')
    Name_system=strcat(system,'_cont');
    % here the same function is called for the continuation system
    % (already including the powers of the parameter)
    interactive_constructor
    
else
    error('this option has not been considered yet')
end

if length(X0)==2*n_nodes+1
    % already Fourier coeffs
    if size(X0,1)==length(X0)
        X0=X0.';
    end
    xBarGen0=Xi_vector([1/period0],X0);
else
    % starting point x(0)
    xBarGen0=solve_system(alpha0,n_nodes,X0,period0);
end

clear coefs_linear

coefs_linear0{3}(1)=0;
for i=1:alpha0.size_scalar
    coefs_linear0{3}(i)=-coefs_linear0{1}(i,:)*xBarGen0.scalar.'-...
        sum(sum(squeeze(coefs_linear0{2}(i,:,:)).*xBarGen0.vector,1));
end

coefs_linear=coefs_linear0;

if bool_plot
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
    if xBarGen0.size_vector==3
        plot3(xBarGen0,'r');
    elseif xBarGen0.size_vector==2
        plot2(xBarGen0,'r');
    else
        plot(xBarGen0,'r');
    end
    hold on
end
temp_intlab=use_intlab;
use_intlab=0;
% NEWTON on the smaller system (parameter fixed)
[x0]=Newton_Xi(xBarGen0,alpha0,coefs_linear0,max_Newton_iter,res_Newton);

% parameter is added
x0.size_scalar=x0.size_scalar+1;
x0.scalar=[x0.scalar,var0];

% non-squared deivative is computed
coefs_linear{1}(:,end+1)=0;
coefs_linear0{1}(:,end+1)=0;
full_DH0=Function_derivative(x0,alpha_coef,coefs_linear0,0);
% flag 0 to have non-square output

% nullspace of the derivative is computed
x_dot0=null(full_DH0);
if bool_invert_dir
    x_dot0=-x_dot0;
end

x_dot0= Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot0,x0.size_scalar,x0.size_vector,...
    x0.nodes)));
x_dot0=x_dot0/norm(x_dot0);

% add one equation for archlength continuation
coefs_linear0{1}=[-x_dot0(1:2).';[0,0]];
coefs_linear0{2}(2,:,:)=coefs_linear{2}(1,:,:);
coefs_linear0{2}(1,:,:)=reshape(-x_dot0(3:end),2*n_nodes+1,x0.size_vector).';
coefs_linear0{3}(2)=coefs_linear{3}(1);
coefs_linear0{3}(1) = Xi_vec2vec(x0).'*x_dot0;

DH0=Function_derivative(x0,alpha_coef,coefs_linear0);


delta_vec=zeros(ITER,1);
Imin_vec=zeros(ITER,1);
Imax_vec=zeros(ITER,1);
xDelta=zeros(ITER,1);
num_nodes_vec=zeros(ITER,1);
Y0_vec=zeros(ITER,x0.size_scalar+x0.size_vector);
Z0_vec=zeros(ITER,x0.size_scalar+x0.size_vector);
Z1_vec=zeros(ITER,x0.size_scalar+x0.size_vector);
Z2_vec=zeros(ITER,x0.size_scalar+x0.size_vector);
normX1_vec=zeros(ITER+1,x0.size_scalar+x0.size_vector);
pho_vec=zeros(ITER+1,1);
normX1_vec(1)=max(norm(x0));
pho_vec(1)=x0.scalar(2); % continuation parameter
% EI=eig(DH0);
% [~,ind]=min(abs(imag(EI)));
% [~,ind2]=min(abs(real(EI)));
% min_eig_imag(1)=imag(EI(ind));
% min_eig_real(1)=real(EI(ind2));

x0=Newton_Xi((x0),alpha_coef,coefs_linear0,max_Newton_iter,res_Newton);

first_iter=0;
use_intlab=temp_intlab;
end


