function bool = Xi_vec2vec_compatible(x_Xi, x_vec)
% function bool = Xi_vec2vec_compatible(x_Xi, x_vec)
% 
% This function tests the compatibility between a vector and a Xi_vector
% This function returns a bool that is true if the two elements are
% compatible and false otherwise. The two elements are compatible if the
% length of the vector is equal to the lenght of Xi_vec2vec(x_Xi), or
% considered otherwise x_Xi.size_scalar+x_Xi.size_vector*(x_Xi.nodes*2+1)
% 
% INPUTS
% x_Xi          Xi_vector
% x_vec         standard vector
%
% OUTPUT
% bool          check if coherent
bool=false;
if  x_Xi.size_scalar+x_Xi.size_vector*(x_Xi.nodes*2+1) == max(size(x_vec)) && min(size(x_vec))==1
    bool=true;
end