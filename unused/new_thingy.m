function [s] = new_thingy (X0,var0,period0,system,ITER,n_nodes,nu_input,...
    delta_step,bool_intlab,bool_display,max_tent,delta_increase,delta_decrease,...
    bool_invert_dir,max_Newton_iter,res_Newton,max_time,bool_plot,bool_talk)
%function [s] = new_thingy (X0,var0,period0,system,ITER,n_nodes,nu_input,...
%     delta_step,bool_intlab,bool_display,max_tent,delta_increase,delta_decrease,...
%     bool_invert_dir,max_Newton_iter,res_Newton,max_time,bool_plot,bool_talk)
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
% bool_talk        if 0 no output,if 1 some output, if 2 more output
%                  (DEFAULT depending on length of the solution vector)
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
global dim_bigger_than_nodes % if the dimension is bigger than 2*nodes+1
global talkative        % if talkative==0 no output
                        % if talkative==1 some output
                        % if talkative==2 more output

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
    Max_comp=max_tent;
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
    if ~exist('bool_talk','var') || isempty(bool_talk) 
        if isempty(talkative)
            if length(X0)*(n_nodes)>300
                talkative=2;
            elseif length(X0)*(n_nodes)>150
                talkative=1;
            else
                talkative=0;
            end
        end
    else
        talkative=bool_talk;
    end
    if talkative>0
        fprintf('Starting the code, time %s\n', datestr(now,13))
    end
    if length(X0)>2*n_nodes+1
        dim_bigger_than_nodes=1;
    else
        dim_bigger_than_nodes=0;
    end
    [x0,DH0,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
        max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
        Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec]=compute_elements(X0,var0,period0,system,ITER,n_nodes,...
        bool_invert_dir,max_Newton_iter,res_Newton,bool_plot);
    
    
else
    if exist('period0','var') && ~isempty(period0)
        if ~isnumeric(period0)
            error('If first input string, second input must be numeric, delta')
        end
        delta_new=period0;
    end
    if exist('system','var') && ~isempty(system)
        if mod(system,1)~=0
            error('If first input string, third input must be integer, number of nodes')
        end
        n_nodes_new=system;
    end
    [x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
        max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
        Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,bool_plot,Max_comp,system,delta_increase,delta_decrease]= ...
        load_elements(X0,var0);%,period0,system,n_nodes_new);
    if x0.size_vector>x0.nodes
        dim_bigger_than_nodes=1;
    else
        dim_bigger_than_nodes=0;
    end
    if exist('delta_new','var')
        delta_step=delta_new;
    end
    temp_intlab=use_intlab;
    if exist('n_nodes_new','var')
        use_intlab=0;
        n_nodes=n_nodes_new;
        x0=reshape_Xi(x0,n_nodes);
        coefs_linear=reshape_coefs(coefs_linear,n_nodes);
        coefs_linear0=reshape_coefs(coefs_linear0,n_nodes);
        
        
        [x0]=Newton_Xi(x0,alpha_coef,coefs_linear0,max_Newton_iter,res_Newton);
        x0=symmetrise(x0);
        
        x_dot0=null(Function_derivative(x0,alpha_coef,coefs_linear,0));
        x_dot0=Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot0,x0.size_scalar,x0.size_vector,...
            x0.nodes)));
        x_dot0=x_dot0/norm(x_dot0);
        
        % add one equation for archlength continuation
        coefs_linear0=coefs_linear;
        coefs_linear0{1}=[-x_dot0(1:2).';[0,0]];
        coefs_linear0{2}(2,:,:)=coefs_linear{2}(1,:,:);
        coefs_linear0{2}(1,:,:)=reshape(-x_dot0(3:end),2*x0.nodes+1,x0.size_vector).';
        coefs_linear0{3}(1) = Xi_vec2vec(x0).'*x_dot0;
        coefs_linear0{3}(2)=0;
        coefs_linear0{3}(2)=real(-coefs_linear0{1}(2,:)*x0.scalar.'-...
            sum(sum(squeeze(coefs_linear0{2}(2,:,:)).*x0.vector,1)));
        
        DH0=Function_derivative(x0,alpha_coef,coefs_linear0);
        use_intlab=temp_intlab;
    end
    n_nodes=x0.nodes;
end


if ~ischar(X0) 
    temp_intlab=use_intlab;
    % before starting, test for the best value of n_nodes
%     possible_nodes=floor(linspace(floor(n_nodes/2),n_nodes*2,10));
%      [best_n_nodes,delta_step]=find_best_n_nodes(possible_nodes,x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
%          max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
%          Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_comp,system,delta_increase,delta_decrease);
    best_n_nodes=n_nodes;
    
    use_intlab=0;
    if best_n_nodes~=n_nodes % updates nodes if necessary
        
        n_nodes=best_n_nodes;
        x0=reshape_Xi(x0,n_nodes);
        coefs_linear=reshape_coefs(coefs_linear,n_nodes);
        coefs_linear0=reshape_coefs(coefs_linear0,n_nodes);
        
        
        [x0]=Newton_Xi(x0,alpha_coef,coefs_linear0,max_Newton_iter,res_Newton);
        x0=symmetrise(x0);
        
        x_dot0=null(Function_derivative(x0,alpha_coef,coefs_linear,0));
        x_dot0=Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot0,x0.size_scalar,x0.size_vector,...
            x0.nodes)));
        x_dot0=x_dot0/norm(x_dot0);
        
        % add one equation for archlength continuation
        coefs_linear0=coefs_linear;
        coefs_linear0{1}=[-x_dot0(1:2).';[0,0]];
        coefs_linear0{2}(2,:,:)=coefs_linear{2}(1,:,:);
        coefs_linear0{2}(1,:,:)=reshape(-x_dot0(3:end),2*x0.nodes+1,x0.size_vector).';
        coefs_linear0{3}(1) = Xi_vec2vec(x0).'*x_dot0;
        coefs_linear0{3}(2)=0;
        coefs_linear0{3}(2)=real(-coefs_linear0{1}(2,:)*x0.scalar.'-...
            sum(sum(squeeze(coefs_linear0{2}(2,:,:)).*x0.vector,1)));
        
        DH0=Function_derivative(x0,alpha_coef,coefs_linear0);
    end
    use_intlab=temp_intlab;
end

tic;
[success, s, x1]=run_algorithm(x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
    max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
    Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_comp,system,delta_increase,delta_decrease);

if bool_plot && success
    %scrsz = get(0,'ScreenSize');
    %figure('Position',[1 scrsz(4)/2 scrsz(3)*3/4 scrsz(4)/2])
    if xBarGen0.size_vector==3
        plot3(x1,'b');
    elseif xBarGen0.size_vector==2
        plot2(x1,'b');
    else
        plot(x1,'b');
    end
    hold on
end

if talkative>0
    fprintf('Code ends, time %s\n', datestr(now,13))
end

return
end











