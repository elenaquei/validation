function coefs_out=reshape_coefs(coefs_in, new_nodes)
% function coefs_out=reshape_coefs(coefs_in, new_nodes)
%
% this function takes as input a 3-cell array, with dimensions as in the
% linear coefficients and reshape the coefficients for a bigger or smaller
% Xi_vector, taking as input the new requested number of nodes

%global use_intlab

coefs_out=coefs_in; %this settle the first and third cell, that will not be changed

old_nodes=(length(coefs_in{2})-1)/2;

if old_nodes<new_nodes
    
    pad=zeros(size(coefs_in{2},1),size(coefs_in{2},2),new_nodes-old_nodes);
    %if ~use_intlab
    coefs_out{2}=cat(3,pad,cat(3,coefs_in{2},pad));
    %else
    %    coefs_out{2}=intvalCAT(3,pad,intvalCAT(3,coefs_in{2},pad));
    %end
elseif old_nodes==new_nodes
    return
else
    diff=old_nodes-new_nodes;
    %if mod(diff,2)~=0
    %    error('requested operation is impossible')
    %end
    coefs_out{2}=coefs_in{2}(:,:,diff+1:end-diff);
end