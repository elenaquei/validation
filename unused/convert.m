function i = convert ( i_iter , deg_f)
% function i = convert ( i_iter , deg_f)

if i_iter > prod(deg_f+1)-1
    warning ('i_iter too big');
end
if i_iter <=0
    error('index must be strictly positive');
end
% i_iter=i_iter -1;
% i = zeros(length(deg_f), 1);
% for n = 1:length(deg_f)
%     remainder = mod(i_iter, deg_f(n));
%     i_iter = floor(i_iter /deg_f(n));
%     i ( n ) = remainder+1;
% end

f = @(k,n,d)mod( floor( k / prod( d(1:n-1))),d(n));
d = deg_f +1;
i = 0*deg_f;
for n =1:length(deg_f)
    i(n) = f(i_iter,n,d);
end