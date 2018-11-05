function [s] = huge_thingy (X0,var0,period0,system,ITER,n_nodes,nu_input,...
    delta_step,bool_intlab,bool_display,max_tent,delta_increase,delta_decrease,...
    bool_invert_dir,max_Newton_iter,res_Newton,max_time,bool_plot)
%function [s] = huge_thingy (X0,var0,period0,system,ITER,n_nodes,nu,...
%    delta_step,bool_intlab,bool_display,max_tent,delta_increase,delta_decrease,...
%    max_Newton_iter,res_Newton,max_time)
% 
% INPUTS
% X0               value of variables at time 0 OR Fourier series for full
%                  [0,period]
% var0             continuation variable
% period0          period of solution (approximate)
% system           either name of the system or coeffs of the system
% ITER             number of iterations (DEFAULT 100)
% Nnodes           number of nodes of the solution (DEFAULT 100)
% nu_input         weight of the l^1_nu space (DEFAULT 11)
% delta_step       length of first step - after adaptative stepsizes 
%                  (DEFAULT 10-3)
% bool_intlab      computations in Intlab or not (DEFAULT 0)
% bool_display     display residual of Newton, all the bounds, Imin and Imax
%                  (DEFAULT 0) 
% max_tent         maximum number of tries for every step (DEFAULT 10)
% delta_increase   if computation succeeds well, stepsize increase by this
%                  value (DEFAULT 1.3)
% delta_decrease   if computation fails, stepsize decreased by this value
%                  (DEFAULT delta_increase)
% bool_invert_dir  invert direction of following the branch (DEFAULT 0)
% max_Newton_iter  maximum number of Newton iterations (DEFAULT 60)
% res_Newton       maximum residual of Newton (DEFAULT 10-14)
% max_time         maximum running time of the code (DEFAULT inf)
% bool_plot        plots first and last orbits (DEFAULT 0)
%
% OR INPUTS
% s_in             string, file where previous results are saved
% ITER             requested number of iterations
% delta_step       stepsize (default: old one)
% n_nodes          number of nodes (default: old one)
%
% OUTPUT
% s                string, file where the results are saved

global nu              % weghted ell^1_\nu norm (FLOAT >=1)
global use_intlab      % if intlab will be used or not (BOOL)
global first_run       % distinguish from the first time this function runs (BOOL)
global Display         % if solutions will be plotted or not (BOOL)

if ~ischar(X0)
    
    first_run=0;
    if ~exist('n_nodes','var')|| isempty(n_nodes)
        n_nodes=100;
    end
    if ~exist('nu_input','var')|| isempty(nu_input)
        nu_input=1.1;
    end
    nu=nu_input;
    if ~exist('delta_step','var')|| isempty(delta_step)
        delta_step=10^-3;
    end
    if ~exist('bool_intlab','var')|| isempty(bool_intlab)
        bool_intlab=1;
    end
    use_intlab=bool_intlab;
    if ~exist('bool_display','var')|| isempty(bool_display)
        bool_display=0;
    end
    Display=bool_display;
    if ~exist('max_tent','var')|| isempty(max_tent)
        max_tent=10;
    end
    Max_Comp=max_tent;
    if ~exist('delta_increase','var')|| isempty(delta_increase)
        delta_increase=1.3;
    end
    if ~exist('delta_decrease','var')|| isempty(delta_decrease)
        delta_decrease=delta_increase;
    end
    if ~exist('bool_invert_dir','var')|| isempty(bool_invert_dir)
        bool_invert_dir=0;
    end
    if ~exist('max_Newton_iter','var')|| isempty(max_Newton_iter)
        max_Newton_iter=60;
    end
    if ~exist('res_Newton','var')|| isempty(res_Newton)
        res_Newton=10^-14;
    end
    if ~exist('max_time','var')|| isempty(max_time)
        max_time=inf;
    end
    if ~exist('bool_plot','var')|| isempty(bool_plot)
        bool_plot=0;
    end
    end_time=0;
    delta=delta_step;
    maxiter=max_Newton_iter;
    min_res=res_Newton;
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
        % here te same function is called for the continuation system
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
    
    previous_iter.Y=[];
    previous_iter.Z1=[];
    
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
    
    good_consecutive_runs=0;

    iter=0;
else
    ITER_new=var0;
    if exist('period0','var') && ~isempty(period0)
        if ~isnumeric(period0)
            error('If first input string, second input must be numeric, delta')
        end
        delta_new=period0;
    end
    if exist('system','var') && ~isempty(system)
        if ~isinteger(system)
            error('If first input string, third input must be integer, number of nodes')
        end
        n_nodes_new=system;
    end
    load(X0);
    ITER=ITER_new+ITER;
    if exist('delta_new','var')
        delta_step=delta_new;
    end
    if exist('n_nodes_new','var')
        x_dot0=Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot0,x0),n_nodes_new));
        x0=reshape_Xi(x0,n_nodes_new);
    end
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
    clear ITER_new
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
end

% % before starting, test for the best value of n_nodes
[best_n_nodes]=find_best_n_nodes(x0,delta_step,coefs_linear,coefs_linear0,alpha_coef,max_Newton_iter,res_Newton);
% use_intlab=0;
% if best_n_nodes~=n_nodes % updates nodes if necessary
%     n_nodes=best_n_nodes;
%     x0=reshape_Xi(x0,n_nodes);
%     coefs_linear=reshape_coefs(coefs_linear,n_nodes);
%     x_dot0=null(Function_derivative(x0,alpha_coef,coefs_linear,0));
%     coefs_linear0=reshape_coefs(coefs_linear0,n_nodes);
%     DH0=Function_derivative(x0,alpha_coef,coefs_linear0);
% end

tic;
while iter<ITER
    iter=iter+1;
    if toc>max_time
        break
    end
    
    % update x1, coefs_linear and x_dot1 without intlab
    use_intlab=0;
    
    [x1,x_dot1,DH1,coefs_linear1] = update(x0,delta_step,x_dot0,...
        alpha_coef,coefs_linear,max_Newton_iter,res_Newton);
    
    use_intlab = temp_intlab;
    
    % validate (x_iter-1 , x_iter)
    [flag,Imin,Imax,previous_iter,Y0,Z0,Z1,Z2,new_step] = radii_polynomials_cont(x0,x1,DH0,DH1,alpha_coef,...
        coefs_linear0, coefs_linear1,previous_iter);
    
    delta_step=delta_step*new_step;
    if new_step>1
        fprintf('Stepsize increased to %1.3d \n',delta_step)
    elseif new_step<1
        fprintf('Stepsize decreased to %1.3d \n',delta_step)
    end
    delta=delta_step;
    revalidation_if_failure;
    delta_step=delta;
    
    comunication_to_user_and_storage;
    
    % if the validation was particularly successfull, the stepsize is incresed
%     if flag==2 && good_consecutive_runs>=5
%         good_consecutive_runs=0;
%         delta_step=delta_step*delta_increase;
%         fprintf('Stepsize increased to %1.3d \n',delta_step)
%     elseif good_consecutive_runs==iter
%         delta_step=delta_step*delta_increase^2;
%         fprintf('Stepsize increased to %1.3d \n',delta_step)
%     end
    
    
    % update of the variables for the next loop
    x0=x1;
    x_dot0=x_dot1;
    DH0=DH1;
    coefs_linear0=coefs_linear1;
    
end

if bool_plot
    %scrsz = get(0,'ScreenSize');
    %figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
    if xBarGen0.size_vector==3
        plot3(xBarGen0,'b');
    elseif xBarGen0.size_vector==2
        plot2(xBarGen0,'b');
    else
        plot(xBarGen0,'b');
    end
    hold on
end


% clear some elements that are useless at this point
clear xBarGen0
clear use_last_run
clear period0
clear coefs_short
clear do_all
clear flag
clear full_DH0
clear i
clear Imax
clear Imin
clear ind
clear ind2
clear nu_input
clear previous_iter
clear sflag
clear X0
clear xDelta
clear xShort0


end_time=end_time+toc;

% s=sprintf('simulations/%s_validation_%d',Name_system,ITER);
s=sprintf('simulations/%s_validation_%d_%d',Name_system,ITER,x0.nodes);
save(s)

