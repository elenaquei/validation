function [F] = F1_simpler_Hopf(X,phi)

lambda=X(1);
beta=X(2);

a1=X(3); b1=X(4); a2=X(5); b2=X(6);
Phi1=[a1;b1]; Phi2=[a2;b2];

x1=0; x2=0;

Df=[[0 1];[-1+2*x1 lambda]];

Phi=Phi1+1i*Phi2;

tempF=Df*Phi-1i*beta*Phi;

F=[phi'*Phi1;phi'*Phi2-1;real(tempF);imag(tempF)];

end

