function [fn,DxFn,Dalphafn] = fn_Rossler(x,alpha)

% from the article
% Hopf Bifurcation in a New Four-Dimensional Hyperchaotic System?
% by 
% LI Xin and YAN Zhen-Ya

A=alpha; 
B = 4;
C = B;
D = 0.04;
E = 1.4;

if D ==0 || B ~= C || B == sqrt(2*E)
    error('chosen parameters wrong')
end

x1 = x(1); y = x(2); z = x(3); w = x(4);


fn=[A*(x1-y)-y*z+w;
    -B*y+x1*z
    -C*z+D*x1+x1*y
    -E*(x1+y)];

DxFn = [A, -A-z, -y, 1
    z, -B, x1, 0
    D+y, x1, -C, 0
    -E, -E, 0,0];
    
Dalphafn= [(x1-y); 0; 0; 0];

end