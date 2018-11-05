function [DFm,DF]=Function_derivative(x,alpha,coefs_linear,flag_square)
% function DFm=Function_derivative(x,alpha,coefs_linear,flag_square)
%
% INPUT:
% x                 Xi_vector, approximate solution;
% alpha             coefs, ODEs equations;
% flag              (DEFAULT=1) the returned matrix will be square
%
% OUTPUT:
% DFm               complex matrix, derivative of the ODE system in x.

global use_intlab

if nargin<3
    error('Too little inputs')
elseif nargin>4
    warning('Too many inputs')
elseif nargin==3
    flag_square=1;
end
if flag_square==0
    %warning('Function_derivative called with flag=0. This flag is for non-square derivatives')
    %flag=1;
end

DF=Xi_matrix(x);
DF.scalar_scalar=coefs_linear{1};
DF.scalar_vector=coefs_linear{2};
K=-x.nodes:x.nodes;
e=eye(DF.size_scalar);

if use_intlab
    DF.scalar_scalar=intval(coefs_linear{1});
    DF.scalar_vector=intval(coefs_linear{2});
    DF.vector_scalar=intval(DF.vector_scalar);
    %DF.vector_vector=intval(DF.vector_vector);
    K=intval(K);
    e=intval(e);
end

for ii=1:DF.size_vector
    for jj=1:DF.size_scalar
        for j=1:alpha.non_zero_el{ii}
            if use_intlab
                temp=alpha.value{ii}(j)*prod(x.scalar.^((alpha.powers_scalars{ii}(:,j)-e(:,jj)).'))...
                    *alpha.powers_scalars{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j));
            else
                temp=alpha.value{ii}(j)*prod(x.scalar.^(alpha.powers_scalars{ii}(:,j).'-e(jj,:)))...
                    *alpha.powers_scalars{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j));
            end
            if j==1
                DF.vector_scalar(ii,:,jj)=DF.vector_scalar(ii,:,jj)+...
                    K.*temp;
            else
                temp=reshape(temp,[1,max(size(temp)),1]);
                DF.vector_scalar(ii,:,jj)=DF.vector_scalar(ii,:,jj)+...
                    temp;
            end
        end
    end
end

e_vec=eye(DF.size_vector);
e_scal=eye(DF.size_scalar);
D=K;%2*pi*x.scalar(1);
size_VV2=size(DF.vector_vector);

Nnodes_new=(size_VV2(3)-1)/2;
if ~use_intlab
    size_vector_scalar=size(DF.vector_scalar);
    size_vector_scalar(2)=Nnodes_new-DF.nodes;
    DF.vector_scalar=cat(2,zeros(size_vector_scalar),...
        DF.vector_scalar,zeros(size_vector_scalar));
    size_scalar_vector=size(DF.scalar_vector);
    size_scalar_vector(3)=Nnodes_new-DF.nodes;
    DF.scalar_vector=cat(3,zeros(size_scalar_vector),...
        DF.scalar_vector,zeros(size_scalar_vector));
else
    DF.vector_vector=intval(DF.vector_vector);
    if Nnodes_new-DF.nodes > 0
        DF.vector_scalar=intvalCAT(2,intvalCAT(2,intval(zeros(DF.size_vector,Nnodes_new-DF.nodes,DF.size_scalar)),...
            DF.vector_scalar),intval(zeros(DF.size_vector,Nnodes_new-DF.nodes,DF.size_scalar)));
        DF.scalar_vector=intvalCAT(3,intvalCAT(3,intval(zeros(DF.size_scalar,DF.size_vector,Nnodes_new-DF.nodes)),...
            DF.scalar_vector),intval(zeros(DF.size_scalar,DF.size_vector,Nnodes_new-DF.nodes)));
    end
end
DF.nodes=Nnodes_new;
Nnodes_all=DF.nodes*2+1;

for ii=1:DF.size_vector % vector function considered
    for jj=1:DF.size_vector % variable considered
        
        for j=1:alpha.non_zero_el{ii}
            
            if j==1 && ii==jj
                continue
            end
            
            temp=alpha.value{ii}(j)*prod(x.scalar.^(alpha.powers_scalars{ii}(:,j).'))...
                *alpha.powers_vectors{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j)-e_vec(:,jj),flag_square);
            if length(temp)<size_VV2(3)
                %if length(temp)~= 0.5*(size_VV2(3)+1)
                %    error('debug')
                %end
                difference=-(length(temp)-size_VV2(3))/2;
                temp=[zeros(1,difference),temp,zeros(1,difference)];
            elseif length(temp)>size_VV2(3)
                mid_temp=(length(temp)+1)/2;
                mid_size=(size_VV2(3)-1)/2;
                temp=temp(mid_temp-mid_size:mid_temp+mid_size);
            end
            if length(temp)~=size_VV2(3)
                error('debug2')
            end
            temp=reshape(temp,[1,1,max(size(temp))]);
            DF.vector_vector(ii,jj,:)=DF.vector_vector(ii,jj,:)+...
                temp;
        end
    end
end

% TO BE TESTED 
%DF2=Function_derivative_tilde(x,alpha,coefs_linear,0);% flag to have small output
%DF3=Function_derivative_tilde2(x,alpha,coefs_linear,0);
%if DF2~=DF3
%    error('here');
%end
DF2=DF;
DF2.vector_scalar=0*DF2.vector_scalar;
for ii=1:DF2.size_vector
    for jj=1:DF2.size_scalar
        for j=1:alpha.non_zero_el{ii}
            if use_intlab
                temp=alpha.value{ii}(j)*prod(x.scalar.^((alpha.powers_scalars{ii}(:,j)-e_scal(:,jj)).'))...
                    *alpha.powers_scalars{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j));
            else
                temp=alpha.value{ii}(j)*prod(x.scalar.^(alpha.powers_scalars{ii}(:,j).'-e_scal(jj,:)))...
                    *alpha.powers_scalars{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j));
            end
            if j==1
                DF2.vector_scalar(ii,:,jj)=DF.vector_scalar(ii,:,jj)+...
                    K.*temp;
            else
                temp=reshape(temp,[1,max(size(temp)),1]);
                DF2.vector_scalar(ii,:,jj)=DF.vector_scalar(ii,:,jj)+...
                    temp;
            end
        end
    end
end

if flag_square % the output is square
    DFm2=Xi_mat2mat(DF2);
    DFm=Xi_mat2mat(DF);
    if any(abs(DFm2-DFm)>eps)
        error('SUPER BIG MISTAKE HERE!') % seems like it's working :)
    end
    size_scalar_col=DF.size_scalar;
    size_scalar_row=DF.size_scalar;
else
    DFm=Xi_mat2mat2(DF,DF.size_scalar-1);
    size_scalar_col=DF.size_scalar;
    size_scalar_row=DF.size_scalar-1;
end

for ii=1:DF.size_vector
    DFm((1+size_scalar_row+Nnodes_all*(ii-1)):(size_scalar_row+Nnodes_all*(ii)),...
        (1+size_scalar_col+Nnodes_all*(ii-1)):(size_scalar_col+Nnodes_all*(ii)))=...
        DFm((1+size_scalar_row+Nnodes_all*(ii-1)):(size_scalar_row+Nnodes_all*(ii)),...
        (1+size_scalar_col+Nnodes_all*(ii-1)):(size_scalar_col+Nnodes_all*(ii)))+...
        alpha.value{ii}(1)*diag(D);
end

return
