% ONE_XI
function  Xi=one_Xi(size_scal,size_vec,nodes)
% Xi   Xi_vector of given size, that has 1-values everywhere
global use_intlab

Xi=Xi_vector(size_scal,size_vec,nodes);
if use_intlab
    Xi.scalar=intval(ones(size(Xi.scalar)));
    Xi.vector=intval(ones(size(Xi.vector)));
else
    Xi.scalar(:)=1;
    Xi.vector(:,:)=1;
end
end
% ONE_XI