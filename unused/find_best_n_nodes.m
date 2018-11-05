function [best_nodes,best_delta]=find_best_n_nodes(possible_nodes,x0,~,delta_step,~,alpha_coef,coefs_linear,coefs_linear0,...
    max_Newton_iter,res_Newton,first_iter,~,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
    Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_comp,system,delta_increase,delta_decrease)
% function [best_nodes,best_delta]=find_best_n_nodes(possible_nodes,x0,~,delta_step,~,alpha_coef,coefs_linear,coefs_linear0,...
%     max_Newton_iter,res_Newton,first_iter,~,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
%     Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_comp,system,delta_increase,delta_decrease)

global use_intlab      % if intlab will be used or not (BOOL)
global talkative

if talkative>0
    disp('Entering find_best_n_nodes')
end
temp_intlab=use_intlab;

%n_nodes=x0.nodes;
%possible_nodes=floor(linspace(floor(n_nodes/2),n_nodes*2,10));

best_nodes=Inf;
best_vel=0;
best_delta=0;
for test_nodes=possible_nodes
    use_intlab=0;
    x0_test=reshape_Xi(x0,test_nodes);
    coefs_linear_test=reshape_coefs(coefs_linear,test_nodes);
    coefs_linear0_test=reshape_coefs(coefs_linear0,test_nodes);
    
    
    [x0_test]=Newton_Xi(x0_test,alpha_coef,coefs_linear0_test,max_Newton_iter,res_Newton);
    x0_test=symmetrise(x0_test);
    
    x_dot0=null(Function_derivative(x0_test,alpha_coef,coefs_linear_test,0));
    x_dot0=Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot0,x0_test.size_scalar,x0_test.size_vector,...
        x0_test.nodes)));
    x_dot0=x_dot0/norm(x_dot0);
    
    % add one equation for archlength continuation
    coefs_linear0_test=coefs_linear_test;
    coefs_linear0_test{1}=[-x_dot0(1:2).';[0,0]];
    coefs_linear0_test{2}(2,:,:)=coefs_linear_test{2}(1,:,:);
    coefs_linear0_test{2}(1,:,:)=reshape(-x_dot0(3:end),2*x0_test.nodes+1,x0_test.size_vector).';
    coefs_linear0_test{3}(1) = Xi_vec2vec(x0_test).'*x_dot0;
    coefs_linear0_test{3}(2)=0; 
    coefs_linear0_test{3}(2)=real(-coefs_linear0_test{1}(2,:)*x0_test.scalar.'-...
        sum(sum(squeeze(coefs_linear0_test{2}(2,:,:)).*x0_test.vector,1)));
    
    DH0=Function_derivative(x0_test,alpha_coef,coefs_linear0_test);
    
    %disp(res);
    use_intlab=temp_intlab;
    
    [flag]=radii_polynomials(x0_test,alpha_coef,coefs_linear0_test);
    
    if flag==0
        continue
    end
    
    use_intlab=0;
    
    tic;
    use_intlab=temp_intlab;
    try
        if talkative>1
            fprintf('\nTesting continuation with %d nodes, try %d out of %d\n\n',test_nodes,find(possible_nodes==test_nodes) , length(possible_nodes));
        end
        [success, ~, ~,new_delta]=run_algorithm(x0_test,DH0,delta_step,x_dot0,alpha_coef,coefs_linear_test,coefs_linear0_test,...
            max_Newton_iter,res_Newton,first_iter,first_iter+1,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
            Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_comp,system,delta_increase,delta_decrease);
        time_run=toc;
        if success && new_delta/time_run >best_vel
            best_nodes = test_nodes;
            best_vel=0.95*new_delta/time_run;
            best_delta=new_delta;
        end
        if success
            figure(101)
            plot(test_nodes, new_delta/time_run,'*')
            hold on;
        end
    catch
        fprintf('failed with %d nodes\n',x0_test.nodes)
    end
end

use_intlab=temp_intlab;

if best_nodes==Inf
    error('problem here')
end
if talkative>0
    disp('Exiting find_best_n_nodes')
end
end