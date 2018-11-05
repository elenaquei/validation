function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
    DalphaalphaxFn,DxxxFnV] = derivatives_hyper(x,a,v)
%function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
%    DalphaalphaxFn,DxxxFnV] = derivatives_hyper(x,a,v)
%

e = 2;
b=1;
c=b;
d=10;

X = x(1); Y = x(2); Z = x(3); W = x(4);

DxFn =[ a    -a-Z   -Y   1
    Z    -b    X   0
    d+Y   X    -c  0
    -e    -e   0   0];

Dalphafn = [X-Y
    0
    0
    0];

if nargout<3
    return
end

v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV = [0      -v3      -v2       0
    v3      0      v1         0
    v2      v1     0          0
    0       0      0          0];

DalphaxFn = zeros(4);

DalphaxxFn=zeros(4);
DalphaalphaFn=zeros(4,1);
DalphaalphaxFn=zeros(4);

DxxxFnV=zeros(4);