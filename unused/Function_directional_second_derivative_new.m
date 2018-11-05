 function DDH=Function_directional_second_derivative_new(xBar,alpha,dir1,dir2)
%  function DDH=Function_directional_second_derivative_new(xBar,alpha,dir1,dir2)
%
% INPUT:
% xBar   Xi_vector
% alpha  full_problem
% dir1   Xi_vector
% dir2   Xi_vector
%
% OUTPUT:
% DDH     Xi_vector

global use_intlab;

DDH=Xi_vector();
DDH.size_scalar=alpha.scalar_equations.num_equations;
DDH.size_vector=alpha.vector_field.n_equations;

if use_intlab
    DDH=intval(DDH);
end


size_scal=xBar.size_scalar;

end_size_vec=xBar.nodes*alpha.vector_field.deg_vector;

% zero padding starts here
halfsize_vec=(xBar.nodes-1)*alpha.vector_field.deg_vector+2;
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
x_conv=get_ifft(xBar); %verifyfft_in(xBar.vector(1,:),-1,1);%
%x_conv=intval(zeros(xBar.size_vector,length(x_conv1)));
%dir1_conv=x_conv;%zeros(xBar.size_vector,length(x_conv1));
%dir2_conv=x_conv;%zeros(xBar.size_vector,length(x_conv1));
%x_conv(1,:)=x_conv1;
dir1_conv=get_ifft(dir1);%verifyfft_in(dir1.vector(1,:),-1,1);
dir2_conv=get_ifft(dir2);%verifyfft_in(dir2.vector(1,:),-1,1);
%for i=2:xBar.size_vector
%    x_conv(i,:)=verifyfft_in(xBar.vector(i,:),-1,1);
%    dir1_conv(i,:)=verifyfft_in(dir1.vector(i,:),-1,1);
%    dir2_conv(i,:)=verifyfft_in(dir2.vector(i,:),-1,1);
%    % taking the long inverse Fourier transform 
%end

%% computation 

%%% second derivative of F(x)
temp_DDF=zeros(alpha.vector_field.n_equations,length(x_conv)); %%% JUST WORKS IF DIM < NODES

if isintval(alpha.vector_field.value) || isintval(x_conv)
    temp_DDF = intval(temp_DDF);
end

for i=1:alpha.vector_field.n_equations % equation
    for j=1:alpha.vector_field.n_terms(i) % element of the equation
        d=[alpha.vector_field.power_scalar{i}(:,j).', alpha.vector_field.power_vector{i}{j}.'];
        const=alpha.vector_field.value{i}(j);
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
                temp_DDF(i,:)= temp_DDF(i,:) + const*d(k)*(d(n)-e(k,n))*...
                    prod(xBar.scalar.^(power_scalar.'))*...
                    product_FFT(x_conv,power_vector).*...
                    temp_dir1.*...
                    temp_dir2;
            end
        end
    end
end


%%% second derivative of G(x)
temp_DDG=zeros(alpha.scalar_equations.num_equations,length(x_conv));
if isintval(temp_DDF) 
    temp_DDG = intval(temp_DDG);
elseif isintval(alpha.scalar_equations.polynomial_equations.value(1))
    temp_DDF = intval(temp_DDF);
    temp_DDG = intval(temp_DDG);
end
    

for i=1:alpha.scalar_equations.number_equations_pol % equation
    % alpha.scalar_equations.number_equations_lin+
    for j=1:alpha.scalar_equations.polynomial_equations.n_terms(i) % element of the equation
        d=[alpha.scalar_equations.polynomial_equations.power_scalar{i}(:,j).', alpha.scalar_equations.polynomial_equations.power_vector{i}{j}.'];
        const=alpha.scalar_equations.polynomial_equations.value{i}(j);
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
                temp_DDG(alpha.scalar_equations.number_equations_lin+i,:)= ...
                    temp_DDG(alpha.scalar_equations.number_equations_lin+i,:) + const*d(k)*(d(n)-e(k,n))*...
                    prod(xBar.scalar.^(power_scalar.'))*...
                    product_FFT(x_conv,power_vector).*...
                    temp_dir1.*...
                    temp_dir2;
            end
        end
    end
end

%% conversion back
temp_DDH = [temp_DDG;temp_DDF];
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


% reshape
DDH=reshape_Xi(DDH,end_size_vec);
