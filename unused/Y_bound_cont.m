function [Yvector,new_Y,Y_s]=Y_bound_cont(A0,A1,xBar0,xBar1,alpha_coef,coefs_linear0,coefs_linear1,old_Y)
% function Yvector=Y_bound_cont(A0,A1,xBar0,xBar1,alpha,coefs_linear0,coefs_linear1)
%
% INPUT
% A0,A1          complex matrices, approximations of the inverse of the 
%                function of the problem;
% xBar0,xBar1    Xi_vector, numerical solutions of the ODE system;
% alpha          coefs, exact data of the problem;
% coefs_linear   coefficients of the linear scalar part.
%
% OUTPUT
% Yvector   a standard vector of length = size_vector+size_scalar, the
%           values of this vector come from the Y bound in the associated pdf (we
%           refer to rigorous numerics for analytic solution of D.E..... for the
%           notation, mainly)
% new_Y     a standard vector of length = size_vector+size_scalar, Y_bound
%           (continuation excluded) useful for next run
% Y_s       Y0+ s^2 YDelta, dependency on the stepsize of Y (rough bound),
%           saved as standard vector of three elements [Y0, 0, YDelta]

global nu
global Display
global use_intlab

if nargin==8
    fullY0=old_Y;
else
    [~,fullY0]=Y_bound(A0,xBar0,alpha_coef,coefs_linear0);
end

[~,fullY1]=Y_bound(A1,xBar1,alpha_coef,coefs_linear1);
new_Y=fullY1;

Y0=cnorm_Xi_vector(max(abs(fullY0),abs(fullY1)),nu);


xBarDelta=xBar0-xBar1;
%
% (!) requires use_intlab anyway (!)
% tic
xBarS=interval_Xi(xBar0,xBar1);
min_coef1=min(coefs_linear1{1},coefs_linear0{1});
max_coef1=min(coefs_linear1{1},coefs_linear0{1});
coefs_linearS{1}=infsup(min_coef1,max_coef1);
min_coef2=min(coefs_linear1{2},coefs_linear0{2});
max_coef2=min(coefs_linear1{2},coefs_linear0{2});
coefs_linearS{2}=infsup(min_coef2,max_coef2);


xBarS= reshape_Xi(xBarS,xBarS.nodes*(alpha_coef.deg_vector-1));
xBarDelta=reshape_Xi(xBarDelta,xBarS.nodes);
if alpha_coef.deg_vector>2
    pad=zeros(size(coefs_linearS{2},1),size(coefs_linearS{2},2),xBar0.nodes*(alpha_coef.deg_vector-2));
    if ~use_intlab
        coefs_linearS{2}=cat(3,pad,cat(3,coefs_linearS{2},pad));
    else
        coefs_linearS{2}=intvalCAT(3,pad,intvalCAT(3,coefs_linearS{2},pad));
    end
end

temp_intlab=use_intlab;
use_intlab=1;
DH_xs2 = Function_derivative(xBarS,alpha_coef,coefs_linearS);


Adiff=extend_approximate_inverse(A1,xBar0.size_scalar,...
    xBar0.size_vector,xBar0.nodes,xBarS.nodes)-...
    extend_approximate_inverse(A0,xBar0.size_scalar,...
    xBar0.size_vector,xBar0.nodes,xBarS.nodes);
Y2=2*cnorm_Xi_vector(vec2Xi_vec(Adiff*DH_xs2*Xi_vec2vec(xBarDelta),...
    xBarDelta.size_scalar,xBarDelta.size_vector,xBarDelta.nodes),nu);

% toc


%     % meaningful change that still doesn't work
%     % should improve the computational time of Y2
%     %%
%     disp('DEBUGGING MODE: maybe faster Y')
%     xBarS=interval_Xi(xBar0,xBar1);
%     xBarDelta=xBar0-xBar1;
% 
%     min_coef1=min(real(coefs_linear1{1}),real(coefs_linear0{1}))+...
%         1i*min(imag(coefs_linear1{1}),imag(coefs_linear0{1}));
%     max_coef1=min(real(coefs_linear1{1}),real(coefs_linear0{1}))+...
%         1i*min(imag(coefs_linear1{1}),imag(coefs_linear0{1}));
%     coefs_linearS{1}=infsup(min_coef1,max_coef1);
% 
%     min_coef2=min(real(coefs_linear1{2}),real(coefs_linear0{2}))+...
%         1i*min(imag(coefs_linear1{2}),imag(coefs_linear0{2}));
%     max_coef2=min(real(coefs_linear1{2}),real(coefs_linear0{2}))+...
%         1i*min(imag(coefs_linear1{2}),imag(coefs_linear0{2}));
%     coefs_linearS{2}=infsup(min_coef2,max_coef2);
% 
%     %nodes_new=xBar1.nodes*alpha_coef.deg_vector;
%     DH_xs3=Function_directional_first_derivative(xBarS,alpha_coef,coefs_linearS,xBarDelta);
%     ADelta_big=reshape_A(A1-A0,xBar1.size_scalar,xBar1.size_vector,xBar1.nodes,...
%         DH_xs3.nodes,xBar1.scalar(1));%-...
%     %reshape_A(,xBar1.size_scalar,xBar1.size_vector,xBar1.nodes,...
%     %DH_xs3.nodes,xBar0.scalar(1));
% 
%     test= Xi_vec2vec(DH_xs3);
% 
%     if any(test ~= DH_xs2*Xi_vec2vec(reshape_Xi(xBarDelta,xBarDelta.nodes*(alpha_coef.deg_vector-1))))
%         error('as predicted, troubles');
%     end
%     Y2=cnorm_Xi_vector(vec2Xi_vec( ADelta_big *Xi_vec2vec(DH_xs3),...
%         xBarDelta.size_scalar,xBarDelta.size_vector,DH_xs3.nodes),nu);
%     %% who knows




if temp_intlab 
    As = infsup (   min(inf(real(A0)),inf(real(A1))), max(sup(real(A0)),sup(real(A1))) )+...
        1i*infsup (   min(inf(imag(A0)),inf(imag(A1))), max(sup(imag(A0)),sup(imag(A1))) );
else
    As = infsup( min(real(A0),real(A1)),max(real(A0),real(A1)) )+...
        1i*infsup( min(imag(A0),imag(A1)),max(imag(A0),imag(A1)) );
end

DDH_xs = Function_directional_second_derivative(xBarS,alpha_coef,xBarDelta,xBarDelta);
% expand A0 to size of vec(DDH_xs)
nodes_new=DDH_xs.nodes;
As=extend_approximate_inverse(As,xBar0.size_scalar,xBar0.size_vector,xBar0.nodes,nodes_new);

Y3=cnorm_Xi_vector(vec2Xi_vec(As*Xi_vec2vec(DDH_xs),...
    xBarDelta.size_scalar,xBarDelta.size_vector,nodes_new),nu);

use_intlab=temp_intlab;

YDelta=(Y2+Y3)/8;
Y_s=[0*Y0, YDelta,0*Y0,Y0];

Yvector=Y0+YDelta;

if any(Yvector>1) && Display
    fprintf('Y_cont computed\n');fprintf(' %d\n',Yvector);
    error('Y_cont is bigger than 1, no interval found')
elseif any(isnan(Yvector))
    error('Y_cont is NaN, big troubles!')
elseif Display
    fprintf('Y_cont computed, %d\n',Yvector);
end

