function y=powers(x,alpha,coefs_linear,flag)
% y=powers(x,alpha,coefs_linear,flag)
% alpha is in the class coefs
% x is in the class Xi_vector
% coefs_linear is the 2-cell having the coefficients of the linear scalar
% part of the system
% the output y is in the class Xi_vector
% flag (DEFAULT=1)   if flag=1, returns a Xi_vector of the same size of x
%                    else returns a Xi_vector of maximum size
global use_intlab

if nargin<4
    flag=1;
end
k=-x.nodes:x.nodes;

if use_intlab
    k=intval(k);
end
% scalar part
if x.size_scalar ~= alpha.size_scalar
    error('dimensions must agree')
end
y=Xi_vector(x);
if isempty(coefs_linear{3}) % temporary correction
    %c=real(sum(x.vector(1,:)));
    for i=1:alpha.size_scalar
        coefs_linear{3}(i)=-coefs_linear{1}(i,:)*x.scalar.'-sum(sum(squeeze(coefs_linear{2}(i,:,:)).*x.vector,1));
    end
    %c=coefs_linear{1}*x.scalar.'+temp;
end

if ~use_intlab
%     temp=(0);
%     for i=1:alpha.size_scalar
%         temp=sum(sum(squeeze((coefs_linear{2}(i,:,:))).*x.vector,1));
%     end
%     y.scalar=(coefs_linear{3})-((coefs_linear{1})*(x.scalar.')).'-temp;
    temp=zeros(size(coefs_linear{3}));%zeros(alpha.size_scalar);
    for i=1:alpha.size_scalar
        temp(i)=sum(sum(squeeze(coefs_linear{2}(i,:,:)).*x.vector,1));
    end
    y.scalar=coefs_linear{3}+(coefs_linear{1}*(x.scalar.')).'+temp;
else
    temp=intval(zeros(size(coefs_linear{3})));
    for i=1:alpha.size_scalar
        temp(i)=sum(sum(squeeze(intval(coefs_linear{2}(i,:,:))).*x.vector,1));
    end
    y.scalar=intval(coefs_linear{3})+(intval(coefs_linear{1})*(x.scalar.')).'+temp;
end


% vector part
if x.size_vector ~= alpha.size_vector
    error('dimensions must agree')
end
if flag==1
    size_vec=size(x.vector,2);
else
    size_vec=(size(x.vector,2)-1)*alpha.deg_vector+1;
end

halfsize_vec=(size_vec-1)/2;
y.vector=zeros(x.size_vector,size_vec);
if use_intlab
    y.vector=intval(y.vector);
end
y.nodes=halfsize_vec;

for i=1:y.size_vector %
    for j=1:alpha.non_zero_el{i}
        if j==1
            temp=k.*alpha.value{i}(j)*prod((x.scalar.').^alpha.powers_scalars{i}(:,j))...
                .*powers_FFT(x.vector,alpha.powers_vectors{i}(:,j),flag);
           
        else
            %y.vector(i,:)=y.vector(i,:)+...
            temp=alpha.value{i}(j)*prod((x.scalar.').^alpha.powers_scalars{i}(:,j))...
                .*powers_FFT(x.vector,alpha.powers_vectors{i}(:,j),flag);
        end
        
        if length(temp)<size_vec
            mid_temp=(length(temp)-1)/2;
            %if length(temp)~= 0.5*(size_vec+1)% NOT VALID HERE
            %    error('debug')
            %end
            if ~use_intlab
                temp=[zeros(1,halfsize_vec-mid_temp),temp,zeros(1,halfsize_vec-mid_temp)];
            else
                
                temp=[intval(zeros(1,halfsize_vec-mid_temp)),temp,intval(zeros(1,halfsize_vec-mid_temp))];
            end
        elseif length(temp)>size_vec
            
            if mod(length(temp),2)
            mid_temp=(length(temp)+1)/2;
            temp=temp(mid_temp-halfsize_vec:mid_temp+halfsize_vec);
            
            else
                
            mid_temp=(length(temp))/2;
            temp=temp(mid_temp-halfsize_vec:mid_temp+halfsize_vec);
            
            end
        end
        %y.vector(i,:)=y.vector(i,:)+tem;
        y.vector(i,:)=y.vector(i,:)+temp;
    end
end
end