function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,DalphaalphaxFn,DxxxFnV] = derivatives_Brusselator(x,alpha,v)

%%% This example is taken from the book
%%%   Kuznetsov, Yuri A. Elements of applied bifurcation theory. 
%%%   Second edition. 112. Springer-Verlag, New York, 1998.
%%%   Exercise (4) (Hopf bifurcation in the Brusselator) of Section 3.6 on
%%%   page 105

A=sqrt(2); B=alpha;
x1=x(1); x2=x(2);

DxFn = [-(B+1)+2*x1*x2      x1^2
        B-2*x1*x2           -x1^2];
    
Dalphafn= [-x1;x1];

if nargout>2

v1=v(1); v2=v(2);

DxxFnV= [ 2*x2*v1 + 2*x1*v2      2*x1*v1
        -2*x2*v1-2*x1*v2        -2*x1*v1];
    
DalphaxFn = [-1,0;1,0];

DalphaxxFn=[0, 0; 0, 0];
DalphaalphaFn=[0;0];
DalphaalphaxFn=zeros(2);

DxxxFnV=[2*v1+2*v2   2*v1;
    -2*v1-2*v2       -2*v1];
end

end

