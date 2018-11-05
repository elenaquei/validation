function x_of_t = Fourrier_of_t(c,t) 

k = (length(c) -1) /2;
if mod(k,1)~=0
    error('length of c must be odd')
end
K= -k:k;
x_of_t = sum(c * exp(1i*2*pi*K*t));

return