function [best_n_nodes]=find_best_n_nodesOLD(x0,delta_step,coefs_linear,coefs_linear0,alpha_coef,max_Newton_iter,res_Newton)
global use_intlab
% ALGORITHM:
% try out for nodes = [1/2, 2/3, 1, 3/2, 2]*n_nodes
% test for the best outcome 
%               BEST: such that delta/time is higher
% if one of the extrema, try outhers as well
% go on with the best
temp_intlab=use_intlab;
previous_iter.Y=[];
previous_iter.Z1=[];
n_nodes=x0.nodes;
nodes_test=floor(1/2*n_nodes):2*n_nodes;
delta_test=0*nodes_test;
time_test=0*nodes_test;
for i=1:length(nodes_test)
    tested_n_nodes=(nodes_test(i));
    
    use_intlab = 0;
    x_test=reshape_Xi(x0,tested_n_nodes);
    
    coefs_linear_test=reshape_coefs(coefs_linear,tested_n_nodes);
    coefs_linear0_test=reshape_coefs(coefs_linear0,tested_n_nodes);
    
    [x_test,~,res]=Newton_Xi(x_test,alpha_coef,coefs_linear0_test,max_Newton_iter,res_Newton);
    x_test=symmetrise(x_test);
    %disp(res);
    use_intlab=temp_intlab;
    
    [flag]=radii_polynomials(x_test,alpha_coef,coefs_linear0_test);
    
    if flag==0
        continue
    end
    
    use_intlab=0;
    x_dot_test=null(Function_derivative(x_test,alpha_coef,coefs_linear_test,0));
    
    x_dot_test=Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot_test,x_test.size_scalar,x_test.size_vector,...
        x_test.nodes)));
    
    coefs_linear0_test{1}=[-x_dot_test(1:2).';[0,0]];
    coefs_linear0_test{2}(1,:,:)=reshape(-x_dot_test(3:end),2*tested_n_nodes+1,x0.size_vector).';
    coefs_linear0_test{3}(1) = Xi_vec2vec(x_test).'*x_dot_test;
    
    DH0_test=Function_derivative(x_test,alpha_coef,coefs_linear0_test);
    
    delta_wigle=2;
    delta_step_test=delta_step/10;
    n_tryies=0;
    while delta_wigle~=1 && n_tryies<10
        tic
        n_tryies=n_tryies+1;
        use_intlab=0;
        delta_step_test=delta_step_test*delta_wigle;
        [x1,~,DH1_test,coefs_linear1_test] = update(x_test,delta_step_test,x_dot_test,...
            alpha_coef,coefs_linear_test,max_Newton_iter,res_Newton);
        
        previous_iter.Y=[];
        previous_iter.Z1=[];
        use_intlab = temp_intlab;
        
        % validate (x_iter-1 , x_iter)
        [flag,~,~,~,~,~,~,~,delta_wigle] = radii_polynomials_cont(x_test,x1,DH0_test,DH1_test,alpha_coef,...
            coefs_linear0_test, coefs_linear1_test,previous_iter);
        time_step_test=toc;
        if flag==0
            time_step_test=Inf;
            delta_wigle=0.5;
        end
    end
    time_test(i)=time_step_test;
    delta_test(i)=delta_step_test;
end
velocity_test=delta_test./time_test;
[~,index_max]=max(velocity_test);
best_n_nodes=nodes_test(index_max)*n_nodes;