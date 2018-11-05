function DF=Function_derivative_tilde(x,alpha,coefs_linear,flag_big)
% function DF=Function_derivative_tilde(x,alpha,coefs_linear,flag_big)
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

global use_intlab

if nargin<3
    error('Too little inputs')
elseif nargin>4
    warning('Too many inputs')
elseif nargin==3
    flag_big=0;
end

% preallocation of the Xi_matrix DF
DF=Xi_matrix(x);
DF.scalar_scalar=coefs_linear{1};
DF.scalar_vector=coefs_linear{2};
e=eye(DF.size_scalar);

if use_intlab
    DF.scalar_scalar=intval(coefs_linear{1});
    DF.scalar_vector=intval(coefs_linear{2});
    DF.vector_scalar=intval(DF.vector_scalar);
    e=intval(e);
end
% part of DF dealing with the scalar derivative of the vector equations
for ii=1:DF.size_vector
    for jj=1:DF.size_scalar
        for j=2:alpha.non_zero_el{ii} % the first element is avoided (tilde)
            if use_intlab
                temp=alpha.value{ii}(j)*prod(x.scalar.^((alpha.powers_scalars{ii}(:,j)-e(:,jj)).'))...
                    *alpha.powers_scalars{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j));
            else
                temp=alpha.value{ii}(j)*prod(x.scalar.^((alpha.powers_scalars{ii}(:,j)-e(:,jj)).'))...
                    *alpha.powers_scalars{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j));
            end
            temp=reshape(temp,[1,max(size(temp)),1]);
            DF.vector_scalar(ii,:,jj)=DF.vector_scalar(ii,:,jj)+...
                temp;
        end
    end
end

e=eye(DF.size_vector);
size_VV2=size(DF.vector_vector);
if flag_big % good resizing of the vector-vector component
    Nnodes_new=(alpha.deg_vector-1)*(size_VV2(3)-1)/2;
    size_VV2(3)=Nnodes_new*2+1;
else
Nnodes_new=(size_VV2(3)-1)/2;
end

% preallocation of the computed ammount of zeros
if ~use_intlab % preallocation 
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
        DF.vector_scalar=intvalCAT(2,intvalCAT(2,...
            intval(zeros(DF.size_vector,Nnodes_new-DF.nodes,DF.size_scalar)),...
            DF.vector_scalar),intval(zeros(DF.size_vector,...
            Nnodes_new-DF.nodes,DF.size_scalar)));
        DF.scalar_vector=intvalCAT(3,intvalCAT(3,...
            intval(zeros(DF.size_scalar,DF.size_vector,Nnodes_new-DF.nodes)),...
            DF.scalar_vector),intval(zeros(DF.size_scalar,DF.size_vector,...
            Nnodes_new-DF.nodes)));
    end
end
DF.vector_vector=zeros(size_VV2);
DF.nodes=Nnodes_new;
if use_intlab
    DF.scalar_vector=intval(DF.scalar_vector);
    DF.vector_scalar=intval(DF.vector_scalar);
    DF.vector_vector=intval(DF.vector_vector);
end

% vector-vector component 
for ii=1:DF.size_vector % vector function considered
    for jj=1:DF.size_vector % variable considered
        
        for j=2:alpha.non_zero_el{ii}
            % derivative computed
            temp=alpha.value{ii}(j)*prod(x.scalar.^(alpha.powers_scalars{ii}(:,j).'))...
                *alpha.powers_vectors{ii}(jj,j).*powers_FFT(x.vector,alpha.powers_vectors{ii}(:,j)-e(:,jj),flag_big);
            
            % good size
            if numel(temp)<size_VV2(3)
                difference=-(length(temp)-size_VV2(3))/2;
                temp=[zeros(1,difference),temp,zeros(1,difference)];
            elseif numel(temp)>size_VV2(3)
                mid_temp=(length(temp)+1)/2;
                mid_size=(size_VV2(3)-1)/2;
                temp=temp(mid_temp-mid_size:mid_temp+mid_size);
            end
            if length(temp)~=size_VV2(3)
                error('debug2')
            end
            % computed element saved
            temp=reshape(temp,[1,1,max(size(temp))]);
            DF.vector_vector(ii,jj,:)=DF.vector_vector(ii,jj,:)+...
                temp;
        end
    end
end

% COMPARISON with Function_derivative_tilde2
% T1=toc;
% tic
DF2=Function_derivative_tilde2(x,alpha,coefs_linear,flag_big);
% T2=toc;
% if T1<T2
%     fprintf('tilde1 faster')
% else
%     fprintf('tilde2 faster')
%     disp(T2-T1);
%     disp(T1)
% end
% if DF2~=DF
%     if any(abs(DF2.vector_vector-DF.vector_vector))>10^-14
%         error('something went wrong again1')
%     end
%     
%     if any(abs(DF2.vector_scalar-DF.vector_scalar))>10^-14
%         error('something went wrong again2')
%     end
% end

return