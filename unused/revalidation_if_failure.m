global use_intlab
global RAD_MAX
global talkative
% if the first try fails, the stepsize will be decreased and the
% computation tempted again. If the computation fails Max_Comp times on
% the same segment, the program crashes
if flag>0
    if RAD_MAX<=Imax
        flag=0;
        fprintf('Problems with RAD_MAX %f',RAD_MAX);
    end
    RAD_MAX=Imax*10;
end
tries=0;
while flag == 0 && tries<Max_Comp
    good_consecutive_runs = 0;
    if iter==1
        % if it's the very first iteration, we change more drastically
        % the stepsize, to ensure convergence
        
        delta=delta/2;
    else
        
        % adding one node
        %         use_intlab=0;
        %         x0=reshape_Xi(x0,x0.nodes+1);
        %         coefs_linear0=reshape_coefs(coefs_linear0,x0.nodes);
        %         coefs_linear=reshape_coefs(coefs_linear,x0.nodes);
        %         x_dot0= Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot0,x0.size_scalar,x0.size_vector,...
        %             x0.nodes-1),x0.nodes));
        %
        %         [x0]=Newton_Xi(x0,alpha_coef,coefs_linear0,maxiter,min_res);
        %         % open discussion: should we keep it or not?
        %
        %         % second question: should we recompute x0dot?
        %
        %         DH0=Function_derivative(x0,alpha_coef,coefs_linear0);
        %use_intlab = temp_intlab;
        previous_iter.Y=[];
        previous_iter.Z1=[];
        
        % decreasing delta a little
        delta=delta/delta_decrease;
    end
    
    % the user is informed of the fact that the stepsize is decreased
    fprintf('The validation did not work, step size is decreased to %1.3d and bounds are recomputed.\n\n', delta)
    
    tries=tries+1;
    % recompute x1, x_dot1 and coefs_linear1 without intlab
    use_intlab=0;
    [x1,x_dot1,DH1,coefs_linear1] = update(x0,delta,x_dot0,...
        alpha_coef,coefs_linear,max_Newton_iter,res_Newton);
    
    use_intlab = temp_intlab;
    
    % validate (x_iter-1 , x_iter)
    [flag,Imin,Imax,previous_iter,Y0,Z0,Z1,Z2,new_step] = radii_polynomials_cont(x0,x1,DH0,DH1,alpha_coef,...
        coefs_linear0, coefs_linear1,previous_iter);
    %     [x1,x_dot1,DH1,coefs_linear1] =update(x0,delta,x_dot0,...
    %         alpha_coef,coefs_linear0,maxiter,min_res);
    %     use_intlab = temp_intlab; %resetting intlab as before
    %
    %     % validate (x_iter-1 , x_iter)
    %     [flag,Imin,Imax,previous_iter,Y0,Z0,Z1,Z2] = radii_polynomials_cont(x0,x1,DH0,DH1,alpha_coef,...
    %         coefs_linear0, coefs_linear1,previous_iter);
end

use_intlab=temp_intlab;
if iter>first_iter+1 && mod(iter,3)==0 %&& 1==0
    n_nodes=x0.nodes;
    possible_nodes=[n_nodes-2,n_nodes,n_nodes+2];
    [best_n_nodes,x1,x_dot1,DH1,Imin,Imax,Y0,Z0,Z1,Z2]=node_testing(...
        x1,x_dot1,coefs_linear1, coefs_linear,alpha_coef, delta_step,...
        max_Newton_iter,res_Newton);
    coefs_linear1=reshape_coefs(coefs_linear,x1.nodes);
    coefs_linear1{1}=[-x_dot1(1:2).';[0,0]];
    coefs_linear1{2}(2,:,:)=coefs_linear1{2}(1,:,:);
    coefs_linear1{2}(1,:,:)=reshape(-x_dot1(3:end),2*x1.nodes+1,x1.size_vector).';
    coefs_linear1{3}(1) = Xi_vec2vec(x1).'*x_dot1;
    coefs_linear1{3}(2)=real(-coefs_linear1{1}(2,:)*x1.scalar.'-...
        sum(sum(squeeze(coefs_linear1{2}(2,:,:)).*x1.vector,1)));
    
    % % % STILL ONE PROBLEM: FIRST SAVE PREVIOUS STEP!
    
    
% %     [best_n_nodes,new_step]=find_best_n_nodes(possible_nodes,x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
% %         max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
% %         Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_Comp,Name_system,delta_increase,delta_decrease);
% % %     %%best_n_nodes=n_nodes-2;
% % %     use_intlab=0;
% % %     if best_n_nodes~=n_nodes % updates nodes if necessary
% % %         previous_iter.Y=[];
% % %         previous_iter.Z1=[];
% % %         if best_n_nodes>n_nodes
% % %             fprintf('\n\n INCREASING NUMBER OF NODES TO %d\n\n',best_n_nodes);
% % %             n_nodes=best_n_nodes;
% % %             x1=reshape_Xi(x1,n_nodes);
% % %             coefs_linear=reshape_coefs(coefs_linear,n_nodes);
% % %             coefs_linear1=reshape_coefs(coefs_linear1,n_nodes);
% % %             
% % %             
% % %             [x1]=Newton_Xi(x1,alpha_coef,coefs_linear1,max_Newton_iter,res_Newton);
% % %             x1=symmetrise(x1);
% % %             
% % %             x_dot1=null(Function_derivative(x1,alpha_coef,coefs_linear,0));
% % %             %             x_dot1=Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot1,x0.size_scalar,x0.size_vector,...
% % %             %                 x1.nodes)));
% % %             %             x_dot1=x_dot1/norm(x_dot1);
% % %             
% % %             angle = atan( imag(x_dot1(2))/real(x_dot1(2)));
% % %             x_dot1 = exp( - 1i * angle) * x_dot1;
% % %             
% % %             x_dot1=Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot1,x1.size_scalar,x1.size_vector,...
% % %                 x1.nodes)));
% % %             x_dot1=x_dot1/norm(x_dot1);
% % %             
% % %             if dot(x_dot1(1:x1.size_scalar),x_dot0(1:x1.size_scalar))<0
% % %                 x_dot1=-x_dot1;
% % %             end
% % %             
% % %             
% % %             % add one equation for archlength continuation
% % %             coefs_linear1=coefs_linear;
% % %             coefs_linear1{1}=[-x_dot1(1:2).';[0,0]];
% % %             coefs_linear1{2}(2,:,:)=coefs_linear{2}(1,:,:);
% % %             coefs_linear1{2}(1,:,:)=reshape(-x_dot1(3:end),2*x1.nodes+1,x1.size_vector).';
% % %             coefs_linear1{3}(1) = Xi_vec2vec(x1).'*x_dot1;
% % %             coefs_linear1{3}(2)=0;
% % %             coefs_linear1{3}(2)=real(-coefs_linear1{1}(2,:)*x1.scalar.'-...
% % %                 sum(sum(squeeze(coefs_linear1{2}(2,:,:)).*x1.vector,1)));
% % %             
% % %             DH1=Function_derivative(x1,alpha_coef,coefs_linear1);
% % %         else % % best_n_nodes < n_nodes
% % %             % one full validation is made, with the idea of decreasing the
% % %             % number of nodes
% % %             % here we'll have two sets of variables, the small and the big. x0
% % %             % starts begin BIG, x1 starts beeing small and padded with zeros to
% % %             % get big
% % %             
% % %             fprintf('\n\n DECREASING NUMBER OF NODES TO %d\n\n',best_n_nodes);
% % %             
% % %             % save previous computation
% % %             comunication_to_user_and_storage;
% % %             
% % %             % start new iteration
% % %             iter=iter+1;
% % %             use_intlab=0;
% % %             
% % %             % BIG dimension
% % %             x0_BIG=x1;
% % %             x_dot0_BIG=x_dot1;
% % %             DH0_BIG=DH1;
% % %             coefs_linear0_BIG=coefs_linear1;
% % %             
% % %             % SMALL dimension (what will be passed on)
% % %             x0_SMALL=reshape_Xi(x0_BIG,best_n_nodes);
% % %             x_dot0_SMALL=Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot1,x0_BIG.size_scalar,x0_BIG.size_vector,...
% % %                 x0_BIG.nodes),x0_SMALL.nodes));
% % %             coefs_linear0_SMALL=reshape_coefs(coefs_linear0_BIG,best_n_nodes);
% % % %             DH0_SMALL=extend_approximate_inverse(DH0_BIG,x0_BIG.size_scalar, x0_BIG.size_vector,x0_BIG.nodes,x0_SMALL.nodes);
% % % %             
% % % %             % compute new step with less nodes
% % %            coefs_linear_SMALL=reshape_coefs(coefs_linear,best_n_nodes);
% % %              [x1_SMALL,x_dot1_SMALL,DH1_SMALL,coefs_linear1_SMALL] = update(x0_SMALL,delta_step/2,x_dot0_SMALL,...
% % %                  alpha_coef,coefs_linear_SMALL,max_Newton_iter,res_Newton);
% % %             
% % %             % padd with zeros
% % %             x1_BIG=reshape_Xi(x1_SMALL,n_nodes);
% % %             coefs_linear1_BIG=reshape_coefs(coefs_linear1_SMALL,n_nodes);
% % %             %             DH1_BIG=extend_numerical_matrix(DH1_SMALL,x1_BIG.size_scalar, x1_BIG.size_vector,x1_SMALL.nodes,x1_BIG.nodes);
% % %             %
% % %             %             x_dot1_BIG=Xi_vec2vec(reshape_Xi(vec2Xi_vec(x_dot1_SMALL,x0_BIG.size_scalar,x0_BIG.size_vector,...
% % %             %                 x0_SMALL.nodes),x0_BIG.nodes));
% % %             %             DH1=Function_derivative(xBar1,alpha_coef,coefs_linear,0);
% % %             %             DH1=[-x_dot1_BIG.';DH1];
% % %             %
% % %             DH1_BIG =Function_derivative(x1_BIG,alpha_coef,coefs_linear1_BIG);
% % %             
% % %             %A1_SMALL=inv(DH1_SMALL);
% % %             %A1_BIG=extend_approximate_inverse(A1_SMALL,x0_BIG.size_scalar, x0_BIG.size_vector,x0_SMALL.nodes,x0_BIG.nodes);
% % %             
% % %             
% % %             use_intlab=temp_intlab;
% % %             
% % %             %if talkative>2
% % %            
% % %             Imin_old=Imin; Imax_old=Imax;
% % %             %end
% % %             previous_iter.Y=[];
% % %             previous_iter.Z1=[];
% % %             DH0_SMALL=Function_derivative(x0_SMALL,alpha_coef,coefs_linear0_SMALL);
% % %             [flag_SMALL,Imin_SMALL,Imax_SMALL,previous_iter_SMALL,Y0_SMALL,...
% % %                Z0_SMALL,Z1_SMALL,Z2_SMALL,new_step_SMALL] = radii_polynomials_cont(x0_SMALL,x1_SMALL,DH0_SMALL,DH1_SMALL,alpha_coef,...
% % %                coefs_linear0_SMALL, coefs_linear1_SMALL,previous_iter);
% % %             [flag,Imin,Imax,previous_iter,Y0,Z0,Z1,Z2,new_step] = radii_polynomials_cont(x0_BIG,x1_BIG,DH0_BIG,DH1_BIG,alpha_coef,...
% % %                 coefs_linear0_BIG, coefs_linear1_BIG,previous_iter);%,[],A1_BIG);
% % %             %if talkative>2
% % %             fprintf('BEFORE DECREASAL [%e, %e]\n',Imin_old,Imax_old)
% % %             fprintf('AFTER DECREASAL [%e, %e]\n',Imin,Imax)
% % %             fprintf('relative difference %e and %e\n\n', (Imin-Imin_old)/Imin, (Imax-Imax_old)/Imax);
% % %             %pause(2)
% % %             %end
% % %             tries=1;
% % %             while flag == 0 && tries<Max_Comp % validation not successful\
% % %                 warning('entered anyway AAAARGGH')
% % %                 fprintf('\n%d\n',iter);
% % %                 good_consecutive_runs = 0; 
% % %                 use_intlab=0; 
% % %                 tries=tries+1;
% % %                 [x1_SMALL,x_dot1_SMALL,DH1_SMALL,coefs_linear1_SMALL] ...
% % %                     = update(x0_SMALL,delta_step/(2*tries),x_dot0_SMALL,...
% % %                     alpha_coef,coefs_linear_SMALL,max_Newton_iter,res_Newton);
% % %                 
% % %                 % padd with zeros
% % %                 x1_BIG=reshape_Xi(x1_SMALL,n_nodes);
% % %                 coefs_linear1_BIG=reshape_coefs(coefs_linear1_SMALL,n_nodes);
% % %                 DH1_BIG =Function_derivative(x1_BIG,alpha_coef,coefs_linear1_BIG);
% % %                 
% % %                 use_intlab=temp_intlab;
% % %                 
% % %                 previous_iter.Y=[];
% % %                 previous_iter.Z1=[];
% % %                 [flag_SMALL] = radii_polynomials_cont(x0_SMALL,x1_SMALL,DH0_SMALL,DH1_SMALL,alpha_coef,...
% % %                     coefs_linear0_SMALL, coefs_linear1_SMALL,previous_iter);
% % %                 [flag,Imin,Imax,previous_iter,Y0,Z0,Z1,Z2,new_step] = radii_polynomials_cont(x0_BIG,x1_BIG,DH0_BIG,DH1_BIG,alpha_coef,...
% % %                     coefs_linear0_BIG, coefs_linear1_BIG,previous_iter);
% % %             end
% % %             
% % %             if flag==0 && tries==Max_Comp
% % %                 % program crashes if the segment is not validated
% % %                 end_time=toc;
% % %                 s=sprintf('%s_crashed_validation_%d',Name_system,iter);
% % %                 save(s)
% % %                 disp(s)
% % %                 error('could not validate the segment')
% % %             end
% % %             % pass the smaller solution, the one with the ideal number of nodes
% % %             
% % %             
% % %             delta_step=delta_step*new_step;
% % %             if new_step>1 && talkative
% % %                 fprintf('Stepsize increased to %1.3d \n',delta_step)
% % %             elseif new_step<1 && talkative
% % %                 fprintf('Stepsize decreased to %1.3d \n',delta_step)
% % %             end
% % %             delta=delta_step;
% % %             
% % %             
% % %             x1_SMALL=x0_SMALL;
% % %             DH1_SMALL =Function_derivative(x0_SMALL,alpha_coef,coefs_linear0_SMALL);
% % %             x_dot1_SMALL=x_dot0_SMALL;
% % %             coefs_linear1_SMALL=coefs_linear0_SMALL;
% % %             
% % %             
% % %             
% % %             x1=x1_SMALL;
% % %             x_dot1=x_dot1_SMALL;
% % %             DH1=DH1_SMALL;
% % %             coefs_linear1=coefs_linear1_SMALL;
% % %             coefs_linear=coefs_linear_SMALL;
% % %             previous_iter.Y=[];
% % %             previous_iter.Z1=[];
% % %         end
% % %     end
    use_intlab=temp_intlab;
end

if tries == 0
    good_consecutive_runs=good_consecutive_runs+1;
else
    good_consecutive_runs = 0;
end