function [Yvector,fullY]=Y_bound(A,xBar,alpha,coefs_linear)
% function Yvector=Y_bound(A,xBar,alpha,coefs_linear)
%
% INPUT
% A      complex matrix, approximation of the inverse of the function of
%        the problem;
% xBar   Xi_vector, numerical solution of the ODE system;
% alpha  coefs, exact data of the problem;
% coefs_linear   coefficients of the linear scalar part.
%
% OUTPUT
% Yvector   a standard vector of length = size_vector+size_scalar, the
%           values of this vector come from the Y bound in the associated pdf (we
%           refer to rigorous numerics for analytic solution od D.E..... for the
%           notation, mainly)

global nu
global use_intlab
global Display 

%Yvector=zeros(xBar.size_vector+xBar.size_scalar,1);

F_bigXi=F_function(xBar,alpha,coefs_linear,0);  % 0 => maximum size
F_bigvec=(F_bigXi);% maximum size         % Xi_vec2vec
F_func_small=F_function(xBar,alpha,coefs_linear);


F_small=Xi_vec2vec(F_func_small);
AF=A*F_small;
AF_Xi=vec2Xi_vec(AF,xBar.size_scalar,xBar.size_vector,xBar.nodes);


onetilde_i_alpha=tilde_Xi(one_Xi(xBar.size_scalar,xBar.size_vector,xBar.nodes*alpha.deg_vector));

if ~use_intlab
    onetilde_i_alpha=sprod_Xi(onetilde_i_alpha,xBar.scalar(1));
else
    onetilde_i_alpha=sprod_Xi(onetilde_i_alpha,xBar.scalar(1));
end

fullY=sum_Xi_vector(AF_Xi,...
    cdivision_Xi_vector(tail(F_bigvec,F_func_small.nodes),onetilde_i_alpha,0),0);
Yvector=cnorm_Xi_vector(fullY,nu);


% stop everything if the result is too big
if any(Yvector>1) && Display
    disp(Yvector);
    error('Y is bigger than 1, no interval found')
elseif Display
    fprintf('Y computed, %d\n',Yvector);
end
