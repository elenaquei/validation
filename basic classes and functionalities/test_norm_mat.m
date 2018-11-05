function [norm_vec]=test_norm_mat(Z_mat,size_scal,size_vec,nodes_col,nodes_row)
%function [norm_vec]=test_norm_mat(Z_mat,size_scal,size_vec,nodes_col,nodes_row)
%
% given size_scal,size_vec, nodes_col, nodes_row
% generalization of the norm of the Xi matrices in case the vertical
% and horizontal number of components is not the same.
global nu 
global use_intlab

if nargin==4
    nodes_row=nodes_col;
elseif nargin<4
    error('not enough input arguments');
end

if use_intlab && ~isintval(Z_mat)
    Z_mat = intval(Z_mat);
end

if ~use_intlab
    norm_scal=sum(abs(Z_mat(1:size_scal,1:size_scal)),2);
else
    norm_scal=sum(sup(abs(Z_mat(1:size_scal,1:size_scal))),2);
end
for i=1:size_vec
    if ~use_intlab
    norm_scal=norm_scal+ sum(abs(Z_mat(1:size_scal,size_scal+(i-1)*(2*nodes_col+1)+(1:2*nodes_col+1)))...
        .*repmat((nu.^ abs(-nodes_col:nodes_col)),[size_scal,1]),2);
    else
        
    norm_scal=norm_scal+ sum(sup(abs(Z_mat(1:size_scal,size_scal+(i-1)*(2*nodes_col+1)+(1:2*nodes_col+1)))...
        .*repmat(nu.^abs(-nodes_col:nodes_col),[size_scal,1])),2);
    end
end
norm_vec=zeros(size_vec,1);
for j=1:size_vec
    for i=1:size_scal
        Zji=Z_mat(size_scal+(j-1)*(2*nodes_row+1)+(1:2*nodes_row+1),i);
        norm_vec(j)=norm_vec(j)+NuNorm(Zji,nu);
    end
    %disp(Z)
    for i=1:size_vec
        Zji=Z_mat(size_scal+(j-1)*(2*nodes_row+1)+(1:2*nodes_row+1),size_scal+(i-1)...
            *(2*nodes_col+1)+(1:2*nodes_col+1));
        norm_vec(j)=norm_vec(j)+InfNuNorm_mat(Zji);
    end
end
norm_vec=[norm_scal;norm_vec];
