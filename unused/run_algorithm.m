function [success,s,x1,delta_step]=run_algorithm(x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
    max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
    Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_Comp,Name_system,delta_increase,delta_decrease)
% function [success,s,x1,delta_step]=run_algorithm(x0,DH0,delta_step,x_dot0,alpha_coef,coefs_linear,coefs_linear0,...
%     max_Newton_iter,res_Newton,first_iter,ITER,delta_vec,Imin_vec,Imax_vec,xDelta,Y0_vec,...
%     Z0_vec,Z1_vec,Z2_vec,num_nodes_vec,normX1_vec,pho_vec,Max_Comp,Name_system,delta_increase,delta_decrease)

global nu              % weghted ell^1_\nu norm (FLOAT >=1)
global use_intlab      % if intlab will be used or not (BOOL)
global first_run       % distinguish from the first time this function runs (BOOL)
global Display         % if solutions will be plotted or not (BOOL)
global talkative

success=0;
temp_intlab=use_intlab;
previous_iter.Y=[];
previous_iter.Z1=[];
good_consecutive_runs=0;
iter=first_iter;

[flag]=radii_polynomials(x0,alpha_coef,coefs_linear0);
if flag==0
    fprintf(datestr(now))
    error('could not validate the segment: first approximation wrong')
end
while iter<ITER
    iter=iter+1;
    %if toc>max_time
    %    break
    %end
    
    % update x1, coefs_linear and x_dot1 without intlab
    use_intlab=0;
    
    [x1,x_dot1,DH1,coefs_linear1] = compute_new_point(x0,delta_step,x_dot0,...
        alpha_coef,coefs_linear,max_Newton_iter,res_Newton);
    
    use_intlab = temp_intlab;
    
    % validate (x_iter-1 , x_iter)
    [flag,Imin,Imax,previous_iter,Y0,Z0,Z1,Z2,new_step] = radii_polynomials_cont(x0,x1,DH0,DH1,alpha_coef,...
        coefs_linear0, coefs_linear1,previous_iter);
    
    delta_step=delta_step*new_step;
    if new_step>1 && talkative
        fprintf('Stepsize increased to %1.3d \n',delta_step)
    elseif new_step<1 && talkative
        fprintf('Stepsize decreased to %1.3d \n',delta_step)
    end
    delta=delta_step;
    revalidation_if_failure;
    delta_step=delta;
    
    if flag==0 && tries==Max_Comp 
        % program crashes if the segment is not validated
        end_time=toc;
        s=sprintf('%s_crashed_validation_%d',Name_system,iter);
        save(s)
        error('could not validate the segment')
    else
        % if the validation succeds, teh user is informed of some of the
        % parameters
        fprintf('%e  + %e i\n', real(x1.scalar(2)), imag(x1.scalar(2)))
        fprintf('norm xDelta %e\n',norm_Xi_vector(x0-x1,nu))
        fprintf('Iteration %d of %d, successful\n', iter, ITER)
        fprintf('Execution time from start: %f sec \n',toc)
        fprintf('Delta: %f\n',delta)
        fprintf('Tries for this step: %d\n\n',tries)
    end
    
    
    
    % if requested, the plot of the nmumerical solution is given 
    if Display 
        plot(x1)
        hold on
    end
    % some elements are stored for further use ( will it be still useful?)
    delta_vec(iter)=delta;
    Imin_vec(iter)=Imin;
    Imax_vec(iter)=Imax;
    xDelta(iter) = norm_Xi_vector(x0-x1,nu);
    Y0_vec(iter,:)=Y0;
    Z0_vec(iter,:)=Z0;
    Z1_vec(iter,:)=Z1;
    Z2_vec(iter,:)=Z2;
    num_nodes_vec(iter)=x0.nodes;
    
    normX1_vec(iter+1,:)=(norm(x1));
    pho_vec(iter+1)=x1.scalar(2);
    
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
success=1;
end_time=toc;

% s=sprintf('simulations/%s_validation_%d',Name_system,ITER);
clear previous_iter
clear DH0
clear DH1
s=sprintf('simulations/%s_validation_%d_%d',Name_system,ITER,x0.nodes);
save(s)
end