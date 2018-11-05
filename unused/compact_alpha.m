function alpha_smaller = compact_alpha(alpha_bigger)
% function alpha_smaller = compact_alpha(alpha_bigger)
%
% INPUT
% alpha_bigger        coefs, possibility of multiple entries with the same
%                     powers
%
% OUTPUT
% alpha_smaller       coefs, no entry with equal powers

alpha_smaller = alpha_bigger;

for j = 1: alpha_bigger.size_vector
    mat = [alpha_bigger.powers_scalars{j} , alpha_bigger.powers_vectors{j}];
    [final , rows, repet ] = unique(mat,'rows');
    alpha_smaller.value{j} = alpha_bigger.value{j}(rows,:);
    u = unique(repet);
    n = histc(repet, u);
    index = find(n>1);
    for k = index
        equal_rows = find(repet==u(k));
        alpha_smaller.value{j}(equal_rows(1)) = sum(alpha_bigger{j}(equal_rows));
    end
    alpha_smaller.powers_scalars{j} = final(:,1:alpha_bigger.size_scalar);
    alpha_smaller.powers_vectors{j} = final(:,alpha_bigger.size_scalar+1:end);
end