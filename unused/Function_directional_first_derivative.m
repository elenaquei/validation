function DH=Function_directional_first_derivative(xBar,alpha,coefs,dir1)
% function DDH=Function_directional_first_derivative(xBar,alpha,dir1)
%
% INPUT:
% xBar   Xi_vector
% alpha  coefs
% dir1   Xi_vector, direction of the derivative
%
% OUTPUT:
% DH     Xi_vector

global use_intlab;

% preallocation
DH=Xi_vector();
DH.size_scalar=xBar.size_scalar;
DH.size_vector=xBar.size_vector;

if use_intlab
    DH=intval(DH);
end


%% scalar side

% if use_intlab
%     DH.scalar = intval(coefs{1})*(dir1.scalar.');
%     for i=1:xBar.size_scalar
%         DH.scalar(i)=DH.scalar(i)+...
%         reshape(squeeze(coefs{2}(i,:,:)).',prod(size(coefs{2}(i,:,:))),1).'*...
%         reshape(dir1.vector.',prod(size((dir1.vector))),1);
%     end
% else
%     
%     DH.scalar = coefs{1}(:,:)*(dir1.scalar.');
%     for i=1:xBar.size_vector
%         DH.scalar=DH.scalar+squeeze(coefs{2}(i,:,:))*dir1.vector;
%     end
% end

if ~use_intlab
    temp=zeros(size(DH.scalar));%zeros(alpha.size_scalar);
    for i=1:alpha.size_scalar
        temp(i)=sum(sum(squeeze(coefs{2}(i,:,:)).*dir1.vector,1));
    end
    DH.scalar=(coefs{1}*(dir1.scalar.')).'+temp;
else
    temp=intval(zeros(size(DH.scalar)));
    for i=1:alpha.size_scalar
        temp(i)=sum(sum(squeeze(intval(coefs{2}(i,:,:))).*dir1.vector,1));
    end
    DH.scalar=(intval(coefs{1})*(dir1.scalar.')).'+temp;
end

%% zero padding

size_scal=xBar.size_scalar;
end_size_vec=xBar.nodes*alpha.deg_vector;

% zero padding starts here
%halfsize_vec=(xBar.nodes-1)*alpha.deg_vector+2;
DH=reshape_Xi(DH,end_size_vec);%halfsize_vec);

if use_intlab
    DH.vector=intval(DH.vector);
end


%% vector part
% first fft

% zero padding
xBar=reshape_Xi(xBar,DH.nodes);
dir1=reshape_Xi(dir1,DH.nodes);

% preallocation for speed
x_conv1=verifyfft_in(xBar.vector(1,:),-1,1);
x_conv=intval(zeros(xBar.size_vector,length(x_conv1)));
dir1_conv=0*x_conv;
x_conv(1,:)=x_conv1;
dir1_conv(1,:)=verifyfft_in(dir1.vector(1,:),-1,1);

% ifft
for i=2:xBar.size_vector
    x_conv(i,:)=verifyfft_in(xBar.vector(i,:),-1,1);
    dir1_conv(i,:)=verifyfft_in(dir1.vector(i,:),-1,1);
    % taking the long inverse Fourier transform 
end

%% computation 

temp_DH=x_conv*0;
for i=1:xBar.size_vector % equation
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
                temp_dir1=prod(dir1.scalar.^(e(k,1:size_scal)));
            else
                temp_dir1=product_FFT(dir1_conv,e(k,size_scal+1:end));
            end
            temp_DH(i,:)= temp_DH(i,:) + const*d(k)*...
                prod(xBar.scalar.^(power_scalar.'))*...
                product_FFT(x_conv,power_vector).*...
                temp_dir1;
            
        end
    end
end

%% conversion back

temp_y_vector=0*temp_DH;
for i=1:xBar.size_vector
    temp_y_vector(i,:)=verifyfft_in(temp_DH(i,:),1,1);
end
nodes_y_vector=(length(temp_y_vector(i,:))-2)/2;
DH.vector=temp_y_vector(:,2:end);
DH.nodes=nodes_y_vector;

% term: 1i omega kx_k 
xBar=xBar.reshape_Xi(nodes_y_vector);
dir1=dir1.reshape_Xi(nodes_y_vector);
k=-nodes_y_vector:nodes_y_vector;
for i=1:xBar.size_vector
    DH.vector(i,:)=DH.vector(i,:)+...
        alpha.value{i}(1)*k.*dir1.vector(i,:);
        %alpha.value{i}(1)*k*dir1.scalar(1).*xBar.vector(i,:)+...
end

% reshape
DH=reshape_Xi(DH,end_size_vec);


