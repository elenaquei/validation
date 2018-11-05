function [DxxFnV] = second_der_hyper(x,v,q)
% function [DxxFnV] = second_der_hyper(x,v,q)
%
e = 2;
b=1;
c=b;
d=10;

X = x(1); Y = x(2); Z = x(3); W = x(4);


v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV = [0      -v3      -v2       0
    v3      0      v1         0
    v2      v1     0          0
    0       0      0          0];
