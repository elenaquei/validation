function DDF = abs_second_derivative_Brusselator(X,~,R)
n=2;
x1=abs(X(1))+R;x2=abs(X(2))+R;
alpha=X(3); beta=X(4); 
v1=abs(X(5:6))+R;v2=abs(X(7:8))+R;

DDF=[2*x2+4*x1-2
    -2*x2-4*x1+2
    2*v1(2)+4*v1(1)+2*(2*x2+4*x1)
    -4*v1(1)-2*v1(2)+2*(-2*x2-4*x1)+4
    2*v2(2)+4*v2(1)+2*(2*x2+4*x1)-4
    -4*v2(1)-2*v2(2)+2*(-2*x1-4*x1)
    0
    0];