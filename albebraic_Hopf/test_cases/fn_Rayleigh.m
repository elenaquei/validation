function [fn] = fn_Rayleigh(x,alpha)

%Rayleigh's equation: x''+ (x')^3 - 2*alpha*x' + x = 0

x1=x(1); x2=x(2); 

fn=[x2;x1+2*alpha*x2-x2^3];

end

