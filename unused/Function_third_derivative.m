function DDDH=Function_third_derivative(xBar,alpha,RAD_MAX)
% function Function_third_derivative(xBar,alpha,RAD_MAX)
%
%
% OUTPUT:
% DDDH     vector

global nu;
global use_intlab;

DDDH=zeros(xBar.size_scalar+xBar.size_vector,1);
%if use_intlab
%    DDDH=intval(DDDH);
% NO: the maximum is wanted anyway
%end

XnormC=cnorm_Xi_vector(xBar,nu);

xANDr=XnormC+RAD_MAX;

size_scal=xBar.size_scalar;

for i=1:xBar.size_vector % equation
    for j=1:alpha.non_zero_el{i} % element of the equation
        d=[alpha.powers_scalars{i}(:,j).', alpha.powers_vectors{i}(:,j).'];
        if use_intlab
            const=sup(abs(alpha.value{i}(j)));
        else
            const=abs(alpha.value{i}(j));
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
                DDDH(i+size_scal)= DDDH(i+size_scal) + abs(const*d(k)*(d(n)-e(k,n))*(d(m)-e(m,k)-e(n,k))*...
                    prod(xANDr.^((d-e(k,:)-e(n,:)-e(m,:)).')));
            end
            end
        end
        
        
    end
end