function DDF = second_derivative_Brusselator(X,~,R)
%
% computation of DDF(y+zR), z\inB_1(0)

n=2;

x1=X(1);x2=X(2);
alpha=X(3); beta=X(4); 
v1=X(5:6);v2=X(7:8);

DDF=[2*x2+4*x1-2
    -2*x2-4*x1+2
    2*v1(2)+4*v1(1)+2*(2*x2+4*x1)
    -4*v1(1)-2*v1(2)+2*(-2*x2-4*x1)+4
    2*v2(2)+4*v2(1)+2*(2*x2+4*x1)-4
    -4*v2(1)-2*v2(2)+2*(-2*x2-4*x1)
    0
    0];

S=sign(DDF);

DDF=abs(DDF+S*R.*[6; 6; 18; 18; 18; 18; 0; 0]);