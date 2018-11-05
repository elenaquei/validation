function c_int = integral_Fourrier_to_F(c,omega)
% function integral = integral_Fourrier_to_F(c,omega)
%
% computuing the integral of the Fourrier series of
% coefficients c and period omega
%  int \sum_k c_k exp(i omega k tau) d tau = \sum_{k~=0} c_k /(i omega
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
if size(K,1)~=size(c,1)
    K=K.';
end


c_tilde= c./ K;
c_tilde(k+1)=0; % taking out the zero-th element

c_int = c_tilde;
c_int (k+1) = - sum( c_tilde );


end