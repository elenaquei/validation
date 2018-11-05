function integral = integral_Fourrier(c,t,omega)
% function integral = integral_Fourrier(c,t,omega)
%
% computuing the integral between 0 and t of the Fourrier series of
% coefficients c and period omega
%  int_0^t \sum_k c_k exp(i omega k tau) d tau = \sum_{k~=0} c_k /(i omega
%  k) exp (i omega k t) - sum_{k~=0} c_k /(i omega k) + c_0 t

k = (length(c) -1) /2;
if mod(k,1)~=0
    error('length of c must be odd')
end
c0 = c(k+1);
K = (1i * omega * [-k:k]);

K(k+1) =0 ; % taking out the zero-th element


if length(K)~=length(c)
    error('something wrong')
end
if size(K,1) ~=size(c,1)
    K=K.';
end
c_tilde= c./ K;
c_tilde(k+1)=0; % taking out the zero-th element

integral = sum( c_tilde.* (exp(K*t)-1) ) + c0 * t;

end