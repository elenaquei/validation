function [fn] = fn_Bautin(x,alpha)

%Bautin's example (p 108 of Kuznetsov's book)
% x1' = x2
% x2' = -x1+ alpha*x2+x1^2+x1*x2+x2^2

x1=x(1); x2=x(2); 
fn=[x2;-x1+ alpha*x2+x1^2+x1*x2+x2^2];

end

