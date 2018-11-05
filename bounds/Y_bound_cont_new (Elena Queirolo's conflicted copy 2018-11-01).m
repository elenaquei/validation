function [Yvector,new_Y,Y_s]=Y_bound_cont_new(A0,A1,xBar0,xBar1,alpha0,alpha1,old_Y)
% function Yvector=Y_bound_cont(A0,A1,xBar0,xBar1,alpha0,alpha1,old_Y)
%
% INPUT
% A0,A1          complex matrices, approximations of the inverse of the
%                function of the problem;
% xBar0,xBar1    Xi_vector, numerical solutions of the ODE system;
% alpha0,alpha1  full_problem
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
global talkative
global use_intlab

if nargin==7
    fullY0=old_Y;
else
    [~,fullY0]=Y_bound_new(A0,xBar0,alpha0);
end

[~,fullY1]=Y_bound_new(A1,xBar1,alpha1);
new_Y=fullY1;

Y0=cnorm_Xi_vector(max(abs(fullY0),abs(fullY1)),nu);


xBarDelta_short=xBar0-xBar1;
%
% (!) requires use_intlab anyway (!)


coefs_linear1 = alpha1.scalar_equations.linear_coef;
coefs_linear0 = alpha0.scalar_equations.linear_coef;
xBarS_short=interval_Xi(xBar0,xBar1);
min_coef1=min(coefs_linear1{1},coefs_linear0{1});
max_coef1=min(coefs_linear1{1},coefs_linear0{1});
coefs_linearS{1}=infsup(min_coef1,max_coef1);
min_coef2=min(coefs_linear1{2},coefs_linear0{2});
max_coef2=min(coefs_linear1{2},coefs_linear0{2});
coefs_linearS{2}=infsup(min_coef2,max_coef2);
min_coef3=min(coefs_linear1{3},coefs_linear0{3});
max_coef3=min(coefs_linear1{3},coefs_linear0{3});
coefs_linearS{3}=infsup(min_coef3,max_coef3);


xBarS= reshape_Xi(xBarS_short,xBarS_short.nodes*(alpha0.vector_field.deg_vector-1));
xBarDelta=reshape_Xi(xBarDelta_short,xBarS.nodes);
if alpha0.vector_field.deg_vector>2
    pad=zeros(size(coefs_linearS{2},1),size(coefs_linearS{2},2),xBar0.nodes*(alpha0.vector_field.deg_vector-2));
    pad = intval(pad);
    coefs_linearS{2}=intvalCAT(3,pad,intvalCAT(3,coefs_linearS{2},pad));
    
end

coefs_linearS{1} = 0*coefs_linearS{1};
coefs_linearS{2} = 0*coefs_linearS{2};
coefs_linearS{3} = 0*coefs_linearS{3};

alphaS = alpha1;
alphaS.scalar_equations.linear_coef = coefs_linearS;

temp_intlab=use_intlab;
use_intlab=1;
%alphaS = reshape(alphaS, xBarS.nodes);
DH_xs2 = derivative_to_matrix(derivative(alphaS,xBarS_short));  % in saddle, here we get NaNs
new_nodes = ((size(DH_xs2,1)-xBar0.size_scalar)/xBar0.size_vector - 1)/2;
xBarDelta_long = reshape_Xi(xBarDelta,new_nodes);

full_DH_xs = DH_xs2*Xi_vec2vec(xBarDelta_long);
result_scalar_eq_Delta = apply(alpha1.scalar_equations,xBarS_short) - apply(alpha0.scalar_equations,xBarS_short);
full_DH_xs(1:alpha1.scalar_equations.num_equations) =result_scalar_eq_Delta;


Adiff=extend_approximate_inverse(A1,xBar0.size_scalar,...
    xBar0.size_vector,xBar0.nodes,new_nodes)-...
    extend_approximate_inverse(A0,xBar0.size_scalar,...
    xBar0.size_vector,xBar0.nodes,new_nodes);
Y2=2*cnorm_Xi_vector(vec2Xi_vec(Adiff*full_DH_xs,...
    xBarDelta_long.size_scalar,xBarDelta_long.size_vector,xBarDelta_long.nodes),nu);



if isintval(A0)
    As = infsup (   min(inf(real(A0)),inf(real(A1))), max(sup(real(A0)),sup(real(A1))) )+...
        1i*infsup (   min(inf(imag(A0)),inf(imag(A1))), max(sup(imag(A0)),sup(imag(A1))) );
else
    As = infsup( min(real(A0),real(A1)),max(real(A0),real(A1)) )+...
        1i*infsup( min(imag(A0),imag(A1)),max(imag(A0),imag(A1)) );
end

DDH_xs = Function_directional_second_derivative(alpha1,xBarS_short,xBarDelta_short,xBarDelta_short);
% addition for the new continuation equation and linear in s scalar eqs
for i = 1:alpha1.scalar_equations.number_equations_lin
    coef_1 = extract_lin_coef(alpha1.scalar_equations, i);
    coef_0 =  extract_lin_coef(alpha0.scalar_equations, i);
    coef_Delta = coef_1-coef_0;
    DDH_xs.scalar(i) = 2*dot(Xi_vec2vec(xBarDelta_short),coef_Delta);
end

% expand A0 to size of vec(DDH_xs)
nodes_new=DDH_xs.nodes;
As=extend_approximate_inverse(As,xBar0.size_scalar,xBar0.size_vector,xBar0.nodes,nodes_new);

Y3=cnorm_Xi_vector(vec2Xi_vec(As*Xi_vec2vec(DDH_xs),...
    xBarDelta.size_scalar,xBarDelta.size_vector,nodes_new),nu);

use_intlab=temp_intlab;

YDelta=(Y2+Y3)/8;
Y_s=[0*Y0, YDelta,0*Y0,Y0];

Yvector=Y0+YDelta;


if any(Yvector>1)
    fprintf('Y_cont computed, %d\n',Yvector);
    error('Y_cont is bigger than 1, no interval found')
elseif any(isnan(Yvector))
    error('Y_cont is NaN, big troubles!')
elseif talkative>1
    fprintf('Y_cont computed, %d\n',Yvector);
elseif talkative>0
    fprintf('Y_cont computed, %d\n',max(Yvector));
end


