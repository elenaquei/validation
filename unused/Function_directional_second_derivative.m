function DDH=Function_directional_second_derivative(xBar,alpha,dir1,dir2)
% function DDH=Function_directional_second_derivative(xBar,alpha,dir1,dir2)
%
% INPUT:
% xBar   Xi_vector
% alpha  coefs
% dir1   Xi_vector
% dir2   Xi_vector
%
% OUTPUT:
% DDH     Xi_vector

global use_intlab;

DDH=Xi_vector();
DDH.size_scalar=xBar.size_scalar;
DDH.size_vector=xBar.size_vector;

if use_intlab
    DDH=intval(DDH);
end

% if isvector(dir1)
% %    dir1=vec2Xi_vec(dir1,xBar.size_scalar,xBar.size_vector,xBar.nodes);
% end
% if isvector(dir2)
% %    dir2=vec2Xi_vec(dir2,xBar.size_scalar,xBar.size_vector,xBar.nodes);
%end

size_scal=xBar.size_scalar;

end_size_vec=xBar.nodes*alpha.deg_vector;

% zero padding starts here
halfsize_vec=(xBar.nodes-1)*alpha.deg_vector+2;
DDH=reshape_Xi(DDH,halfsize_vec);

if use_intlab
    DDH.vector=intval(DDH.vector);
end


%% vector part
% computing

% way to go:
% ifft(xvector)
xBar=reshape_Xi(xBar,DDH.nodes);
dir1=reshape_Xi(dir1,DDH.nodes);
dir2=reshape_Xi(dir2,DDH.nodes);


% preallocation for speed
x_conv1=verifyfft_in(xBar.vector(1,:),-1,1);
x_conv=intval(zeros(xBar.size_vector,length(x_conv1)));
dir1_conv=x_conv;%zeros(xBar.size_vector,length(x_conv1));
dir2_conv=x_conv;%zeros(xBar.size_vector,length(x_conv1));
x_conv(1,:)=x_conv1;
dir1_conv(1,:)=verifyfft_in(dir1.vector(1,:),-1,1);
dir2_conv(1,:)=verifyfft_in(dir2.vector(1,:),-1,1);
for i=2:xBar.size_vector
    x_conv(i,:)=verifyfft_in(xBar.vector(i,:),-1,1);
    dir1_conv(i,:)=verifyfft_in(dir1.vector(i,:),-1,1);
    dir2_conv(i,:)=verifyfft_in(dir2.vector(i,:),-1,1);
    % taking the long inverse Fourier transform 
end

%% computation 

temp_DDH=x_conv*0;
for i=1:xBar.size_vector % equation
    for j=2:alpha.non_zero_el{i} % element of the equation
        d=[alpha.powers_scalars{i}(:,j).', alpha.powers_vectors{i}(:,j).'];
        const=alpha.value{i}(j);
        N=length(d);
        e=eye(N);
        
        % DDH_i = sum_{k=1,d_k>0}^N sum_{j=1,(d-e_k)_n>0}^N
        % (d-e_k)_n d_k x^{d-e_k-e_n} * dir1^e_k * dir2^e_n
        
        for k=1:N
            %if d(k)<=0 
            %    continue
            %end
            for n=1:N
                %if d(n)-e(n,k)<=0
                %    continue
                %end
                power=((d-e(k,:)-e(n,:)).');
                power_scalar=power(1:size_scal);
                power_vector=power(size_scal+1:end);
                if any(power<0)
                    continue
                end
                if k<=size_scal
                    temp_dir1=dir1.scalar(k);%prod(dir1.scalar.^(e(k,1:size_scal)));
                else
                    temp_dir1=dir1_conv(k-size_scal,:);%product_FFT(dir1_conv,e(k,size_scal+1:end));
                end
                if n<=size_scal
                    temp_dir2=dir2.scalar(n);%prod(dir2.scalar.^(e(n,1:size_scal)));
                else
                    temp_dir2=dir2_conv(n-size_scal,:);%product_FFT(dir2_conv,e(n,size_scal+1:end));
                end
                temp_DDH(i,:)= temp_DDH(i,:) + const*d(k)*(d(n)-e(k,n))*...
                    prod(xBar.scalar.^(power_scalar.'))*...
                    product_FFT(x_conv,power_vector).*...
                    temp_dir1.*...
                    temp_dir2;
            end
        end
    end
end

%% conversion back

temp_y_vector=0*temp_DDH;
for i=1:xBar.size_vector
    if all(temp_DDH(i,:)==0)
        continue
    end
    temp_y_vector(i,:)=verifyfft_in(temp_DDH(i,:),1,1);
end
nodes_y_vector=(length(temp_y_vector(i,:))-2)/2;
DDH.vector=temp_y_vector(:,2:end);
DDH.nodes=nodes_y_vector;

% xBar=xBar.reshape_Xi(nodes_y_vector);
% dir1=dir1.reshape_Xi(nodes_y_vector);
% dir2=dir2.reshape_Xi(nodes_y_vector);
% k=-nodes_y_vector:nodes_y_vector;
% for i=1:xBar.size_vector
%     DDH.vector(i,:)=DDH.vector(i,:)+...
%         alpha.value{i}(1)*k*dir1.scalar(1).*dir2.vector(i,:)+...
%         alpha.value{i}(1)*k*dir2.scalar(1).*dir1.vector(i,:);
% end

% reshape
DDH=reshape_Xi(DDH,end_size_vec);
