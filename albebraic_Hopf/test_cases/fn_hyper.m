function [fn,DxFn,Dalphafn] = fn_hyper(x,a)
% function [fn,DxFn,Dalphafn] = fn_hyper(x,a)
%

e = 2;
b=1;
c=b;
d=10;

X = x(1); Y = x(2); Z = x(3); W = x(4);

fn = [ a*X-a*Y-Y*Z+W
    -b*Y+X*Z
    -c*Z+d*X+X*Y
    -e*X-e*Y];

if nargout<2 
    return
end

DxFn =[ a    -a-Z   -Y   1
    Z    -b    X   0
    d+Y   X    -c  0
    -e    -e   0   0];

Dalphafn = [X-Y
    0
    0
    0];