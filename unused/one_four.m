function y= one_four(x)

k = (length(x) -1) /2;
if mod(k,1)~=0
    error('length of c must be odd')
end
y=0*x;
y(k+1)=1;

return