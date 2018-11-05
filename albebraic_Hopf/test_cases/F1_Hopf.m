function [F] = F1_Hopf(X,phi)

lambda=X(1);
beta=X(2);
x1=X(3); x2=X(4);
a1=X(5); b1=X(6); a2=X(7); b2=X(8);
Phi1=[a1;b1]; Phi2=[a2;b2];

f=[x2;-x1*(1-x1)+lambda*x2];
Df=[[0 1];[-1+2*x1 lambda]];

Phi=Phi1+1i*Phi2;

tempF=Df*Phi-1i*beta*Phi;

F=[phi'*Phi1;phi'*Phi2-1;f;real(tempF);imag(tempF)];

end

