function DDDH=Function_third_derivative_new(xBar,alpha,RAD_MAX)
% function Function_third_derivative_Ew(xBar,alpha,RAD_MAX)
%
% INPUTS
% xBar     Xi_vector, numerical solution
% alpha    full_problem, vector field and scalar equations
% RAD_MAX  double, estimation of the maximum validation radius
%
% OUTPUT:
% DDDH     vector

global nu;
global use_intlab;
size_alpha = alpha.scalar_equations.num_equations + alpha.vector_field.n_equations;
DDDH=zeros(size_alpha,1);

XnormC=cnorm_Xi_vector(xBar,nu);

xANDr=XnormC+RAD_MAX;

%% DDG
for i=1:alpha.scalar_equations.number_equations_pol
    %1+alpha.scalar_equations.number_equations_lin : alpha.vector_field.n_equations% equation
    for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
        d=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
        if use_intlab
            const=sup(abs(alpha.scalar_equations.polynomial_equations.value{i}(j)));
        else
            const=abs(alpha.scalar_equations.polynomial_equations.value{i}(j));
        end
        N=length(d);
        e=eye(N);
        
        % DDH_i = sum_{k=1,d_k>0}^N sum_{n=1,(d-e_k)_n>0}^N sum_{m=1,(d-e_k-e_n)_m>0}^N
        % (d-e_k-e_n)_m (d-e_k)_n d_k x^{d-e_k-e_n}
        
        for k=1:N
            if d(k)<0 
                continue
            end
            for n=1:N
                if d(n)-e(n,k)<=0
                    continue
                end
                
            for m=1:N
                if d(m)-e(m,k)-e(n,k)<=0
                    continue
                end
                DDDH(i+alpha.scalar_equations.number_equations_lin)= DDDH(i+alpha.scalar_equations.number_equations_lin)...
                    + abs(const*d(k)*(d(n)-e(k,n))*(d(m)-e(m,k)-e(n,k))*...
                    prod(xANDr.^((d-e(k,:)-e(n,:)-e(m,:)).')));
            end
            end
        end
        
        
    end
end

%% DDF
for i=1 : alpha.vector_field.n_equations% equation
    for j=1:alpha.vector_field.n_terms(i) % element of the equation
        d=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
        if use_intlab
            const=sup(abs(alpha.vector_field.value{i}(j)));
        else
            const=abs(alpha.vector_field.value{i}(j));
        end
        N=length(d);
        e=eye(N);
        
        % DDH_i = sum_{k=1,d_k>0}^N sum_{n=1,(d-e_k)_n>0}^N sum_{m=1,(d-e_k-e_n)_m>0}^N
        % (d-e_k-e_n)_m (d-e_k)_n d_k x^{d-e_k-e_n}
        
        for k=1:N
            if d(k)<0 
                continue
            end
            for n=1:N
                if d(n)-e(n,k)<=0
                    continue
                end
                
            for m=1:N
                if d(m)-e(m,k)-e(n,k)<=0
                    continue
                end
                DDDH(i+alpha.scalar_equations.num_equations)= DDDH(i+alpha.scalar_equations.num_equations)...
                    + abs(const*d(k)*(d(n)-e(k,n))*(d(m)-e(m,k)-e(n,k))*...
                    prod(xANDr.^((d-e(k,:)-e(n,:)-e(m,:)).')));
            end
            end
        end
        
        
    end
end