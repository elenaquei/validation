function DDH=Function_second_derivative(xBar,alpha,RAD_MAX)
% function Function_second_derivative(xBar,alpha,RAD_MAX)  GOOD
%
%
% OUTPUT:
% DDH     vector

global nu;
global use_intlab;

DDH=zeros(xBar.size_scalar+xBar.size_vector,1);
if use_intlab
%    DDH=intval(DDH);
% the output should be the maximum of the interval anyway
end

XnormC=cnorm_Xi_vector(xBar,nu);

xANDr=XnormC+RAD_MAX;

size_scal=xBar.size_scalar;

for i=1:xBar.size_vector % equation
    for j=2:alpha.non_zero_el{i} % element of the equation
        d=[alpha.powers_scalars{i}(:,j).', alpha.powers_vectors{i}(:,j).'];
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
                DDH(i+size_scal)= DDH(i+size_scal) + abs(const*d(k)*(d(n)-e(k,n))*...
                    prod(xANDr.^((d-e(k,:)-e(n,:)).')));
            end
        end
        
    end
end