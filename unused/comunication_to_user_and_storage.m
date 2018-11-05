
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
    
%     EI=eig(DH1);
%     [~,ind]=min(abs(EI));
%     %[~,ind2]=min(abs(real(EI)));
%     min_eig(iter+1)=min(abs(EI));
%     min_eig_imag(iter+1)=imag(EI(ind));
%     min_eig_real(iter+1)=real(EI(ind));
    
%     if pho_vec(iter+1) >0.224799 && delta>10^-6
%        delta=10^-6;
%     end