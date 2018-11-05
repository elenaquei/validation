function [fn] = fn1(x,lambda)

fn=[x(2);-x(1)*(1-x(1))+lambda*x(2)];

end

