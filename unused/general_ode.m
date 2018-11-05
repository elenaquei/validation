function y=general_ode(alpha,t,x)
%function general_ode(alpha,t,x)
%
% this function takes as input the coefficients of the Fourier system and
% the standard time and space coordinates and compute the standard system
%
% INPUT
% alpha    coefs
% t        time (remark: the system is considered to be autonomous)
% x        space
%
% OUTPUT
% F(alpha,t,x)
%
% N.B.: It is considered that for each function the first term stored in
% alpha is the derivative, and therefore it will be ignored.
flag2=0;
xold=x;
% % flag=0;
% % xold=x;
% % if min(size(x))==1
% % else
% %     error('problem with input');
% % 
% %     flag=1;
% %     size_x=max(size(xold))/alpha.size_vector;
% %     x=zeros(2,size_x);
% %     for i=1:alpha.size_vector
% %         x(i,:)=xold(1+(i-1)*size_x:i*size_x);
% %     end
% % end
% % if max(size(x))~=alpha.size_vector
% %     error('dimensions wrong')
% % elseif size(x,1)~=2||size(x,1)~=3
% %     %flag2=1;
% %     %x=x';
% % end
if min(size(x))==1
    flag=1;
    size_x=max(size(xold))/alpha.size_vector;
    x=zeros(2,size_x);
    for i=1:alpha.size_vector
        x(i,:)=xold(1+(i-1)*size_x:i*size_x);
    end
elseif min(size(x))~=alpha.size_vector
    error('dimensions wrong')
elseif size(x,1)~=2
    flag2=1;
    x=x';
end
%%%%%% I don't get it anymore...
y=0*x;

for i=1:alpha.size_vector
    for j=2:alpha.non_zero_el{i}
        temp=1+0*x(i,:);
        for k=1:alpha.size_vector
            temp=temp*x(k,:).^alpha.powers_vectors{i}(k,j);
        end
        y(i,:)=y(i,:)+alpha.value{i}(j)*temp;
    end
end

if flag==1
    y=reshape(y,numel(y),1);
end
if flag2==1
    y=y.';
end
if size(y)~=size(xold)
    if size(y')==size(xold)
        y=y.';
    else
        error('something wrong');
        % they are
        %disp('size are ok')
    end
end