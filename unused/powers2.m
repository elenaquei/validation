function y=powers2(x,alpha,coefs_linear,flag)
% y=powers2(x,alpha,coefs_linear,flag)   %%working properly
% alpha is in the class coefs
% x is in the class Xi_vector
% coefs_linear is the 2-cell having the coefficients of the linear scalar
% part of the system
% the output y is in the class Xi_vector
% flag (DEFAULT=1)   if flag=1, returns a Xi_vector of the same size of x
%                    else returns a Xi_vector of maximum size
global c
global use_intlab

if nargin<4
    flag=1;
end

if x.size_scalar ~= alpha.size_scalar
    error('dimensions must agree')
end
if x.size_vector ~= alpha.size_vector
    error('dimensions must agree')
end

if x.nodes ~=(size(x.vector,2)-1)/2
    error('Input does not agree (x.nodes not mathing size of x.vector)')
end

y=Xi_vector();
y.size_scalar=x.size_scalar;
y.size_vector=x.size_vector;

%% scalar part

%if isempty(c) % temporary correction Not compatible with cont
    %c=real(sum(x.vector(1,:)));
    temp=0;
    for i=1:alpha.size_scalar
        temp=sum(sum(squeeze(coefs_linear{2}(i,:,:)).*x.vector,1));
    end
    c=coefs_linear{1}*x.scalar.'+temp;
%end

if ~use_intlab
    temp=0;
    for i=1:alpha.size_scalar
        temp=sum(sum(squeeze(coefs_linear{2}(i,:,:)).*x.vector,1));
    end
    y.scalar=c-coefs_linear{1}*(x.scalar.')-temp;
else
    temp=intval(0);
    for i=1:alpha.size_scalar
        temp=sum(sum(squeeze(intval(coefs_linear{2}(i,:,:))).*x.vector,1));
    end
    y.scalar=intval(c)-intval(coefs_linear{1})*(x.scalar.')-temp;
end

tic;
%% vector part
% setting constants

if flag==1 % storing the final requested length of the output vector 
    end_size_vec=x.nodes;
else
    end_size_vec=x.nodes*alpha.deg_vector;
end

% zero padding starts here
halfsize_vec=(x.nodes-1)*alpha.deg_vector+2;
y=reshape_Xi(y,halfsize_vec);

if use_intlab
    y.vector=intval(y.vector);
end


%% vector part
% computing

% way to go:
% ifft(xvector)
x=reshape_Xi(x,y.nodes);


% preallocation for speed
x_conv1=verifyfft_in(x.vector(1,:),-1,1);
x_conv=zeros(x.size_vector,length(x_conv1));
if use_intlab
    x_conv=intval(x_conv);
end
x_conv(1,:)=x_conv1;
for i=2:x.size_vector
    x_conv(i,:)=verifyfft_in(x.vector(i,:),-1,1);
    % taking the long inverse Fourier transform 
end

nodes_x_conv=(length(x_conv(i,:))-2)/2;
k=-(nodes_x_conv+1):nodes_x_conv;
if use_intlab
    k=intval(k);
end

% applying function
temp_y=zeros(y.size_vector,length(k));
if use_intlab
    temp_y=intval(temp_y);
end
for i=1:y.size_vector 
    for j=2:alpha.non_zero_el{i}
        temp=alpha.value{i}(j)*prod((x.scalar.').^alpha.powers_scalars{i}(:,j))...
            .*product_FFT(x_conv,alpha.powers_vectors{i}(:,j));
        
        temp_y(i,:)=temp_y(i,:)+temp;
    end
end

%% conclusion
% conversion back

temp_y_vector=0*temp_y;
for i=1:x.size_vector
    temp_y_vector(i,:)=verifyfft_in(temp_y(i,:),1,1);
    if sum(alpha.powers_vectors{i}(:,1))>0
        temp_y_vector(i,:)=temp_y_vector(i,:)+...
            k.*(verifyfft_in(alpha.value{i}(1)*prod((x.scalar.').^alpha.powers_scalars{i}(:,1))...
            .*product_FFT(x_conv,alpha.powers_vectors{i}(:,1)),1,1).');
    else
        temp_y_vector(i,-k(1)+1)=temp_y_vector(i,-k(1)+1)+1;
    end
end
nodes_y_vector=(length(temp_y_vector(i,:))-2)/2;
y.vector=temp_y_vector(:,2:end);
y.nodes=nodes_y_vector;

% reshape
y=reshape_Xi(y,end_size_vec);

end



%%%%%%%%%%%
%           Helping functions
%%%%%%%%%%%





