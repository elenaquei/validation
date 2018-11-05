% CONTINUATION_MAIN
global use_intlab 

% if use_last_run is 1, the elements in the workspace will be used, if it
% is 0 they will be reloaded or reasked for
use_last_run=1;

if ~use_last_run || ~exist('iter','var')
    load_necessary_elements;
end

max_Newton_iter=maxiter;
res_Newton= min_res;


%for iter=1:ITER
while iter<ITER
    iter=iter+1;
    if toc>max_time
        break
    end
    
    if iter==iter_stop || iter==170
        pause(1);
    end
    
    % update x1, coefs_linear and x_dot1 without intlab
    use_intlab=0;
    
    [x1,x_dot1,DH1,coefs_linear1] = update(x0,delta,x_dot0,...
        alpha_coef,coefs_linear,maxiter,min_res);
    %[x1test,x_dot1test,DH1test,coefs_linear1test] = update_2(x0,delta,x_dot0,...
    %    alpha_coef,alpha_short,coefs_short,coefs_linear,maxiter,min_res);

    
    use_intlab = temp_intlab;
    
    % validate (x_iter-1 , x_iter)
    [flag,Imin,Imax,previous_iter,Y0,Z0,Z1,Z2] = radii_polynomials_cont(x0,x1,DH0,DH1,alpha_coef,...
        coefs_linear0, coefs_linear1,previous_iter);
    % flag=1;Imin=0;Imax=0;previous_iter=0;Y0=0;Z0=0;Z1=0;Z2=0;
    revalidation_if_failure;
    
    comunication_to_user_and_storage;
    
    % if the validation was particularly successfull, the stepsize is incresed
    if flag==2 && good_consecutive_runs>=5
        good_consecutive_runs=0;
        delta=delta*delta_increase;
        fprintf('Stepsize increased to %1.3d \n',delta)
    elseif good_consecutive_runs==iter
        delta=delta*delta_increase^2;
        fprintf('Stepsize increased to %1.3d \n',delta)
    end
    
    
    % update of the variables for the next loop
    x0=x1;
    x_dot0=x_dot1;
    DH0=DH1;
    coefs_linear0=coefs_linear1;
    
end

end_time=toc;

s=sprintf('simulations/%s_validation_%d',Name_system,ITER);
s=sprintf('simulations/%s_validation_%d_%d',Name_system,ITER,x0.nodes);
save(s)

return