function DDH=Function_second_derivative_new(xBar,alpha,RAD_MAX)
% function Function_second_derivative_new(xBar,alpha,RAD_MAX)  GOOD
%
% INPUT:
% xBar     Xi_vector
% alpha    full_problem
% RAD_MAX  positive scalar
%
% OUTPUT:
% DDH     vector

global nu;
global use_intlab;


if ~isa(alpha, 'full_problem')
    error('input not suitable')
end



DDH=zeros(alpha.scalar_equations.num_equations+alpha.vector_field.n_equations,1);
if use_intlab
%    DDH=intval(DDH);
% the output should be the maximum of the interval anyway
end

XnormC=cnorm_Xi_vector(xBar,nu);

xANDr=XnormC+RAD_MAX;

size_scal=xBar.size_scalar;

%% here for the polynomial scalar equations

% the derivative is null for hte linear part,
% alpha.scalar_equation.number_equations_lin, but non-null for the scalar
% polynomial part

n_scalar_equations = alpha.scalar_equations.number_equations_lin;

beta = alpha.scalar_equations.polynomial_equations;

for i=1:alpha.scalar_equations.number_equations_pol % equation
    for j=1:beta.n_terms(i) % element of the equation
        if any(size(beta.power_vector{i}{j},2)>1)
            error('Second derivative not yet ready for this')
        end
        d=[beta.power_scalar{i}(:,j).', beta.power_vector{i}{j}.'];
        if use_intlab
            const=sup(abs(beta.value{i}(j)));
        else
            const=abs(beta.value{i}(j));
        end
        N=length(d);
        e=eye(N);
        
        % DDH_i = sum_{k=1,d_k>0}^N sum_{j=1,(d-e_k)_n>0}^N
        % (d-e_k)_n d_k x^{d-e_k-e_n}
        
        for k=1:N
            if d(k)<=0 
                continue
            end
            for n=1:N
                if d(n)-e(n,k)<=0
                    continue
                end
                DDH(i+n_scalar_equations)= DDH(i+n_scalar_equations) + abs(const*d(k)*(d(n)-e(k,n))*...
                    prod(xANDr.^((d-e(k,:)-e(n,:)).')));
            end
        end
        
    end
end

%% here for the vector field part

n_scalar_equations = alpha.scalar_equations.num_equations;
alpha = alpha.vector_field;

for i=1:xBar.size_vector % equation
    for j=1:alpha.n_terms(i) % element of the equation
        if any(size(alpha.power_vector{i}{j},2)>1)
            error('Second derivative not yet ready for this')
        end
        d=[alpha.power_scalar{i}(:,j).', alpha.power_vector{i}{j}.'];
        if use_intlab
            const=sup(abs(alpha.value{i}(j)));
        else
            const=abs(alpha.value{i}(j));
        end
        N=length(d);
        e=eye(N);
        
        % DDH_i = sum_{k=1,d_k>0}^N sum_{j=1,(d-e_k)_n>0}^N
        % (d-e_k)_n d_k x^{d-e_k-e_n}
        
        for k=1:N
            if d(k)<=0 
                continue
            end
            for n=1:N
                if d(n)-e(n,k)<=0
                    continue
                end
                DDH(i+n_scalar_equations)= DDH(i+n_scalar_equations) + abs(const*d(k)*(d(n)-e(k,n))*...
                    prod(xANDr.^((d-e(k,:)-e(n,:)).')));
            end
        end
        
    end
end