 function mat=Xi_mat2mat(Xi_mat)
% function mat=Xi_mat2mat(Xi_mat)
%
% this function takes as input a Xi_matrix and returns a normal matrix,
% such that Xi_mat applied to any Xi_vector is equal to mat*vec, where vec
% is the convertion of the Xi_vector into a vector (as done in Xi_vec2vec).
%

global use_intlab

k=Xi_mat.size_scalar+Xi_mat.size_vector*(2*Xi_mat.nodes+1);

%mat=zeros(k);

mat1=Xi_mat.scalar_scalar;

M=Xi_mat.size_scalar;
N=Xi_mat.size_vector;
F=2*Xi_mat.nodes+1;

mat2=reshape(permute(Xi_mat.scalar_vector,[1,3,2]),M,N*F);

mat3=reshape(permute(Xi_mat.vector_scalar,[2,1,3]),N*F,M);

if ~use_intlab
mat4=zeros(N*F,N*F);
else
mat4=intval(zeros(N*F,N*F));
end
% Xi_mat.vec_vec is actually a NxF vector, such that 
% Xi_mat.vec applyied to Xi_v is the convolution of the two N-vectors. Now we
% want to built a matrix such that the matrix product is equivalent to the
% convolution in all cases

for i=1:N
    for j=1:N
        mat4((i-1)*F+1:i*F,(j-1)*F+1:j*F)=conv2mat(Xi_mat.vector_vector(i,j,:),F);
    end
end

mat=[mat1,mat2;mat3,mat4];

% debugging check
if size(mat,1)~=size(mat,2)
    error('1')
elseif size(mat,1)~=k
    error('2')
end

end
