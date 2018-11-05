function [F] = F_VdP(X,phi)

%Van der Pol's oscillator: x'' - ( alpha - x^2 )*x' + x = 0

alpha=X(1); % Parameter 
beta=X(2); % distance of the complex conjugate eigs from the real axis
x1=X(3); x2=X(4); % the variables from the model
a1=X(5); b1=X(6); a2=X(7); b2=X(8); % 
Phi1=[a1;b1]; Phi2=[a2;b2];

f=[x2;(alpha-x1^2)*x2-x1];
Df=[[0 1];[-1-2*x1*x2 alpha-x1^2]];

Phi=Phi1+1i*Phi2;

tempF=Df*Phi-1i*beta*Phi;

F=[phi'*Phi1;phi'*Phi2-1;f;real(tempF);imag(tempF)];

end

