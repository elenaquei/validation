function [Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector)
% function [Imin,Imax]=find_negative(Z2vector,Z1vector,Z0vector,Yvector)
% 
% this function computes the interval in which the radii polynomial is
% negative
%       P(r) = Yvector + (Z0vector + Z1vector -1) r + Z2vector r^2
%
% INPUT
% Z2vector              real vector, second order term
% Z1vector, Z0vector    real vector, first order term
% Yvector               real vector, zeroth order term
%                 all vectors must be of the same length
% OUTPUT
%  Imin, Imax           lower and upper bound of the interval in which the
%                       radii polynomial

global use_intlab
global RAD_MAX

if use_intlab
    a= intval(Z2vector);
    b = intval(Z1vector)+intval(Z0vector);
    b = b - intval(1);
    c = intval(Yvector);
else
    a = Z2vector;
    b = Z1vector + Z0vector -1;
    c = Yvector;
end

Delta = b.^2 - 4*a.*c;

if any((Delta)<0)
    error('No interval found');
end

Imin = zeros(length(a),1);
Imax = zeros(length(a),1);

for i = 1:length(a)
    if a(i)>0
        Imin(i) = ( - b(i) - sqrt(Delta(i)))./(2*a(i));
        Imax(i) = ( - b(i) + sqrt(Delta(i)))./(2*a(i));
    else
        Imin(i) = -c(i)/b(i);
        Imax(i) = 10;
    end
end

if use_intlab 
    Imin = sup(Imin);
    Imax = inf(Imax);
else
    Imin = max(Imin,0);
end
Imin = max(Imin);
Imax = min(Imax);
if Imin >Imax
    error('No interval found')
elseif Imin>RAD_MAX
    error('Minimum radius bigger than expected radius, used for Z2 bound')
end
Imax = min(Imax, RAD_MAX);

return
end

