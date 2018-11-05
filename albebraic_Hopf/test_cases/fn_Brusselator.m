function [fn,DxFn,Dalphafn] = fn_Brusselator(x,alpha)

%%% This example is taken from the book
%%%   Kuznetsov, Yuri A. Elements of applied bifurcation theory. 
%%%   Second edition. 112. Springer-Verlag, New York, 1998.
%%%   Exercise (4) (Hopf bifurcation in the Brusselator) of Section 3.6 on
%%%   page 105

A=sqrt(2); B=alpha;
x1=x(1); x2=x(2);

fn=[A-(B+1)*x1+x1^2*x2;B*x1-x1^2*x2];

DxFn = [-(B+1)+2*x1*x2      x1^2
        B-2*x1*x2           -x1^2];
    
Dalphafn= [-x1;x1];

end

