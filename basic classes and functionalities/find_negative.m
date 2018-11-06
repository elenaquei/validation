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
end

return
end


% sol=zeros(length(Yvector),2);
% % matrix of coefficients
% Koefs=[Z2vector,Z1vector+Z0vector-1,Yvector];
% 
% 
% %for i=1:length(Yvector)
% %    sol1(i,:)=roots(Koefs(i,:));
% %end
% 
% % computation of solution
% % explicit formula can be used, since the radii polynomail is always second
% % order
% if use_intlab
%     a=intval(Koefs(:,1));b=intval(Koefs(:,2));c=intval(Koefs(:,3));
%     sol(:,1)=inf( (-b-(b.^2-4*a.*c).^0.5) ./ (2*a));
%     sol(:,2)=inf( (-b+(b.^2-4*a.*c).^0.5) ./ (2*a));
% else
%     a=Koefs(:,1);b=Koefs(:,2);c=Koefs(:,3);
%     sol(:,1)= (-b-(b.^2-4*a.*c).^0.5) ./ (2*a);
%     sol(:,2)= (-b+(b.^2-4*a.*c).^0.5) ./ (2*a);
% end
% 
% 
% % 
% index = find(Z2vector~=0);
% if use_intlab
%     a=intval(Koefs(index,1));b=intval(Koefs(index,2));c=intval(Koefs(index,3));
%     sol(index,1)=inf( (-b-(b.^2-4*a.*c).^0.5) ./ (2*a));
%     sol(index,2)=inf( (-b+(b.^2-4*a.*c).^0.5) ./ (2*a));
% else
%     a=(Koefs(index,1));b=(Koefs(index,2));c=(Koefs(index,3));
%     sol(index,1)= (-b-(b.^2-4*a.*c).^0.5) ./ (2*a);
%     sol(index,2)= (-b+(b.^2-4*a.*c).^0.5) ./ (2*a);
% end
% index = find(Z2vector==0);
% if use_intlab
%     a=intval(Koefs(index,1));b=intval(Koefs(index,2));c=intval(Koefs(index,3));
%     sol(index,1)=inf(-c./b);
% else
%     a=(Koefs(index,1));b=(Koefs(index,2));c=(Koefs(index,3));
%     sol(index,1)= -c./b;
% end
% sol(index,2) =Inf;
% 
% 
% 
% % error if any complex interval is found
% if any(imag(sol))
%     Imin=[];Imax=[];
%     error('No interval found')
% end
% 
% % minima and maximum are checked to be positive  // COULD BE BETTER
% min_int=zeros(length(Yvector),1);
% max_int=zeros(length(Yvector),1);
% for i=1:length(Yvector)
%     try
%         min_int(i)=min(sol(i,sol(i,:)>=0));
%     catch
%         error('No convergence interval');
%     end
%     max_int(i)=max(sol(i,:));
% end
% 
% if find(max_int<=0)
%     Imin=[];Imax=[];
%     error('No interval found')
% end
% % general minimum and maximum computed
% Imin=max(min_int);
% Imax=min(max_int);
% 
% if any(size(Imin)==0)
%     error('No interval found')
% end