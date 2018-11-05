function DDDF = third_der_Brusselator(x, v1,v2,v3)
% function DDDF = third_der_Brusselator(x, v1,v2,v3)

a=v1(1); b=v1(2);
c=v2(1); d=v2(2);
e=v3(1); f=v3(2);
x1=x(1); x2=x(2);

DDDF=[2*b*c*e+2*b*d*e+2*a*c*f
    -2*x2*b*c*e-2*a*d*e-2*x2*b*d*e-2*a*c*f-2*x1*b*c*f-2*x1*b*d*f];
