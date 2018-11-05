function [fn] = fn2(x,alpha)

%%% This example is taken from the book
%%%   Kuznetsov, Yuri A. Elements of applied bifurcation theory. 
%%%   Second edition. 112. Springer-Verlag, New York, 1998.
%%%   Example 3.1 (Hopf bifurcation in a predator-prey model) of Section 3.5 on page 101 

r=1; c=2; d=1;
x1=x(1); x2=x(2);

fn=[r*x1*(1-x1)-c*x1*x2/(alpha+x1);-d*x2+c*x1*x2/(alpha+x1)];

end

