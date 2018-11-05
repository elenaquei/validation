function mat=conv2mat(vector,size_mat)
% function mat=conv2mat(vector,size_mat)
%
% INPUT 
% vector    vector
% size_mat  integer scalar, size of output
% OUTPUT
% mat       matrix
%
% takes a input a vector and returns the square matrix created such that
% mat * vec2 = conv(vector,vec2)

m=(length(vector)-1)/2;

vector=reshape(vector,1,length(vector));

% vector padded with zeros if necessary
while size_mat>m
    vector=[0*vector, vector, 0*vector];
    m=(length(vector)-1)/2;
end

% important part selected
row1=vector(m+1:-1:m+2-size_mat);
column1=vector(m+1:m+size_mat);

% toeplitz built-in function used
mat= toeplitz(column1,row1);

% debugging
if size(mat,1)~=size_mat
    error('debugging error');
end

end
