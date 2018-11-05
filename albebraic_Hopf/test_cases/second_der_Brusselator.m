function [DxxFnV] = second_der_Brusselator(x,q1,q2)

%%% This example is taken from the book
%%%   Kuznetsov, Yuri A. Elements of applied bifurcation theory. 
%%%   Second edition. 112. Springer-Verlag, New York, 1998.
%%%   Exercise (4) (Hopf bifurcation in the Brusselator) of Section 3.6 on
%%%   page 105


x1=x(1); x2=x(2);

v1=q1(1); v2=q1(2);

DxxFnV= [ 2*x2*v1 + 2*x1*v2      2*x1*v1
        -2*x2*v1-2*x1*v2        -2*x1*v1]*q2;

end

