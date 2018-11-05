function [DxFn,Dalphafn,DxxFnV,DalphaxFn,DalphaxxFn,DalphaalphaFn,...
    DalphaalphaxFn,DxxxFnV] = derivatives_Rossler(x,alpha,v)

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

DxFn = [A, -A-z, -y, 1
    z, -B, x1, 0
    D+y, x1, -C, 0
    -E, -E, 0,0];
    
Dalphafn= [(x1-y); 0; 0; 0];

if nargout>2

v1=v(1); v2=v(2); v3=v(3); v4 = v(4);

DxxFnV= [0, -v3, -v2, 0
    v3, 0, v1, 0
    v2, v1, 0, 0
    0, 0, 0, 0];
    
DalphaxFn =[1, -1, 0, 0
    0, 0, 0, 0
    0, 0, 0, 0
    0, 0, 0,0];

DalphaxxFn=zeros(4);
DalphaalphaFn=zeros(4,1);
DalphaalphaxFn=zeros(4);

DxxxFnV=zeros(4);
end

end

