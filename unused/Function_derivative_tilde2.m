function DH=Function_derivative_tilde2(x,alpha,coefs,flag_big)
% function DF=Function_derivative_tilde2(x,alpha,coefs_linear,flag_big)
% without the 2 *pi*1i*1tilde term
%
% INPUT
%
% x                     Xi_vector
% alpha                 coeffs
% flag_big=0 (DEFAULT)  the vector_vector component will have just $nodes$ elements
% flag_big=1            the vector_vector component will have $(deg-1)nodes$ elements
%
% OUTPUT
%
% DF                    Xi_matrix

global use_intlab;


if nargin<3
    error('Too little inputs')
elseif nargin>4
    warning('Too many inputs')
elseif nargin==3
    flag_big=0;
end

DH=Xi_matrix(x);

if use_intlab
    DH=intval(DH);
end


%% scalar side

DH.scalar_scalar=coefs{1};
DH.scalar_vector=coefs{2};

%% zero padding

size_scal=x.size_scalar;
if flag_big
    end_size_vec=x.nodes*(alpha.deg_vector-1);
else
    end_size_vec=x.nodes;
end

% zero padding starts here
Nnodes_new=end_size_vec;
DH.vector_vector=zeros(x.size_vector,x.size_vector,Nnodes_new*2+1);
if ~use_intlab % preallocation 
    size_vector_scalar=size(DH.vector_scalar);
    size_vector_scalar(2)=Nnodes_new-DH.nodes;
    DH.vector_scalar=cat(2,zeros(size_vector_scalar),...
        DH.vector_scalar,zeros(size_vector_scalar));
    size_scalar_vector=size(DH.scalar_vector);
    size_scalar_vector(3)=Nnodes_new-DH.nodes;
    DH.scalar_vector=cat(3,zeros(size_scalar_vector),...
        DH.scalar_vector,zeros(size_scalar_vector));
else
    DH.vector_vector=intval(DH.vector_vector);
    DH.vector_scalar=intval(DH.vector_scalar);
    if Nnodes_new-DH.nodes > 0
        DH.vector_scalar=intvalCAT(2,intvalCAT(2,intval(...
            zeros(DH.size_vector,Nnodes_new-DH.nodes,DH.size_scalar)),...
            DH.vector_scalar),intval(zeros(DH.size_vector,...
            Nnodes_new-DH.nodes,DH.size_scalar)));
        DH.scalar_vector=intvalCAT(3,intvalCAT(3,intval(...
            zeros(DH.size_scalar,DH.size_vector,Nnodes_new-DH.nodes)),...
            DH.scalar_vector),intval(zeros(DH.size_scalar,...
            DH.size_vector,Nnodes_new-DH.nodes)));
    end
end
DH.nodes=Nnodes_new;


%% vector part
% first fft

% zero padding
x=reshape_Xi(x,DH.nodes);

% preallocation for speed
x_conv1=verifyfft_in(x.vector(1,:),-1,1);
if use_intlab
    x_conv=intval(zeros(x.size_vector,length(x_conv1)));
end
x_conv(1,:)=x_conv1;
% ifft
for i=2:x.size_vector
    x_conv(i,:)=verifyfft_in(x.vector(i,:),-1,1);
    % taking the long inverse Fourier transform 
end

%% computation 

temp_DH_vv=zeros(x.size_vector,x.size_vector,length(x_conv));
temp_DH_vs=zeros(x.size_vector,length(x_conv),x.size_scalar);%x_conv*0;
if use_intlab
    temp_DH_vv=intval(temp_DH_vv);
    temp_DH_vs=intval(temp_DH_vs);
end
for i=1:x.size_vector % equation
    for j=2:alpha.non_zero_el{i} % element of the equation
        d=[alpha.powers_scalars{i}(:,j).', alpha.powers_vectors{i}(:,j).'];
        const=alpha.value{i}(j);
        N=length(d);
        e=eye(N);
        
        % DH_i = sum_{k=1,d_k>0}^N d_k x^{d-e_k} * dir1^e_k
        
        for k=1:N
            if d(k)<=0
                continue
            end
            
            power=((d-e(k,:)).');
            power_scalar=power(1:size_scal);
            power_vector=power(size_scal+1:end);
            if k<=size_scal
                
            temp_DH_vs(i,:,k)= temp_DH_vs(i,:,k) + const*d(k)*...
                prod(x.scalar.^(power_scalar.'))*...
                product_FFT(x_conv,power_vector);
            %temp_dir1=prod(dir1.scalar.^(e(k,1:size_scal)));
            else
                
            temp_DH_vv(i,k-size_scal,:)= temp_DH_vv(i,k-size_scal,:) + permute(const*d(k)*...
                prod(x.scalar.^(power_scalar.'))*...
                product_FFT(x_conv,power_vector),[3,1,2]);
            end
        end
    end
end

%% conversion back

temp_y_vectorvector=0*temp_DH_vv;
temp_y_vectorscalar=0*temp_DH_vs;
for i=1:x.size_vector
    for k=1:x.size_vector
        temp_y_vectorvector(i,k,:)=verifyfft_in(temp_DH_vv(i,k,:),1,1);
    end
    for k=1:size(temp_DH_vs,3)%x.size_scalar
        temp_y_vectorscalar(i,:,k)=verifyfft_in(temp_DH_vs(i,:,k),1,1);
    end
end
nodes_y_vector=(length(temp_y_vectorvector(i,i,:))-2)/2;
%DH.nodes=nodes_y_vector;

% reshape
nodes_temp_y=(size(temp_y_vectorvector,3)-2)/2;
%check dimensions
if size(temp_y_vectorvector,3)~=size(temp_y_vectorscalar,2)
    error('some big mistake');
end
diff_nodes=nodes_temp_y-end_size_vec;
if diff_nodes<0
    error('something else went wrong');
elseif diff_nodes==0
    DH.vector_vector=temp_y_vectorvector(:,:,2:end);
    DH.vector_scalar=temp_y_vectorscalar(:,2:end,:);
    return
end

% get back to requested size
DH.vector_vector=temp_y_vectorvector(:,:,(2+diff_nodes):(end-diff_nodes));
DH.vector_scalar=temp_y_vectorscalar(:,(2+diff_nodes):(end-diff_nodes),:);

if size(DH.vector_vector,3)~=(end_size_vec*2+1)
    error('debugging mistake');
end
DH.nodes=end_size_vec;
% DH=reshape_Xi(DH,end_size_vec);


