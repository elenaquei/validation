% set of variables loaded and first solution computed
if ~exist('x0','var') || use_last_run==0
    set_parameters;
    
    temp_intlab=use_intlab;
    use_intlab=0;
    
    fprintf('main_Xi_cont\n\n')
    fprintf('\n\nSolving system %s \n\n\n',Name_system);
    
    [x0, x_dot0, ~,xShort0,alpha_coef,coefs_linear,alpha_short,coefs_short] = ...
        compute_solution_and_derivative(Name_system,flag_plot,n_nodes,maxiter,min_res);
end

%if x_dot0(2) >0 % continuation variable increases
    % go in the other direction (sometimes necessary)
    % x_dot0=-x_dot0;
%end

% add to the coefficients of the scalar equations the extra equation 
%     E_s(x) = (x_s - x) \dot \dot{x}_s

if ~exist('coefs_linear0{1}','var') || use_last_run==0
    % add one equation for archlength continuation
    coefs_linear0=coefs_linear;
    
    coefs_linear0{1}=[-x_dot0(1:2).';[0,0]];
    coefs_linear0{2}(2,:,:)=coefs_linear0{2}(1,:,:);
    coefs_linear0{2}(1,:,:)=reshape(-x_dot0(3:end),2*n_nodes+1,x0.size_vector).';
    coefs_linear0{3}(2)=coefs_linear{3}(1);
    coefs_linear0{3}(1) = Xi_vec2vec(x0).'*x_dot0;
end

% squared derivative computed
if ~exist('DH0','var') || use_last_run==0
    DH0=Function_derivative(x0,alpha_coef,coefs_linear0);
end

% clear previous_iter;
if ~exist('previous_iter','var') || use_last_run==0
    previous_iter.Y=[];
    previous_iter.Z1=[];
end
%time count starts here
tic
iter=0;

% preallocation of some vectors to store interesting values during the
% iterations
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


good_consecutive_runs=0;