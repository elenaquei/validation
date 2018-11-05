function [F] = F2_Hopf(X,phi)

%%% This example is taken from the book
%%%   Kuznetsov, Yuri A. Elements of applied bifurcation theory. 
%%%   Second edition. 112. Springer-Verlag, New York, 1998.
%%%   Example 3.1 (Hopf bifurcation in a predator-prey model) of Section 3.5 on page 101 

alpha=X(1);
beta=X(2);

r=1; c=2; d=1;

x1=X(3); x2=X(4);
a1=X(5); b1=X(6); a2=X(7); b2=X(8);
Phi1=[a1;b1]; Phi2=[a2;b2];

f=[r*x1*(1-x1)-c*x1*x2/(alpha+x1);-d*x2+c*x1*x2/(alpha+x1)];

Df11=r*(1-x1)-r*x1-c*x2/(alpha+x1)+c*x1*x2/(alpha+x1)^2;
Df12=-c*x1/(alpha+x1);
Df21=c*x2/(alpha+x1)-c*x1*x2/(alpha+x1)^2;
Df22=-d+c*x1/(alpha+x1);
Df=[[Df11 Df12];[Df21 Df22]];

Phi=Phi1+1i*Phi2;

tempF=Df*Phi-1i*beta*Phi;

F=[phi'*Phi1;phi'*Phi2-1;f;real(tempF);imag(tempF)];

end

