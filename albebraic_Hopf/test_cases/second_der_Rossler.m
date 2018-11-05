function [DxxFnV] = second_der_Rossler(x,v,q2)
% from the article
% Hopf Bifurcation in a New Four-Dimensional Hyperchaotic System?
% by 
% LI Xin and YAN Zhen-Ya

x1=x(1); y=x(2);
z = x(3); w = x(4);

v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV= [0, -v3, -v2, 0
    v3, 0, v1, 0
    v2, v1, 0, 0
    0, 0, 0, 0]*q2;

end