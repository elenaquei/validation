function [Z1vector,new_Z1,Z1_s]=Z1_bound_cont_new(A0,A1,xBar0,xBar1,alpha0,alpha1,Adagger_delta,old_Z1)
% function [Z1vector,new_Z1,Z1_s]=Z1_bound_cont_new(A0,A1,xBar0,xBar1,alpha0,alpha1,Adagger_delta,old_Z1)
% 
% if given, old_Z1 needs to have the following structure:
%   old_Z1.mat (matrix)
%   old_Z1.nodes_col (positive integer)
%   old_Z1.nodes_row (positive integer)
%
% Z1_s        Z1+ s^2 Z1Delta, dependency on the stepsize of Z1 (rough bound),
%             saved as standard vector of three elements [0,Z1Delta,0, Z1]


global use_intlab
global nu
global talkative

N=xBar0.size_vector;
M=xBar0.size_scalar;


if nargin==9 && ~isempty(old_Z1)
    Z_mat0=old_Z1.mat;
    nodes_col0=old_Z1.nodes_col;
    nodes_row0=old_Z1.nodes_row;
else
    [~,Z_mat0,nodes_col0,nodes_row0]=Z1_bound_new(A0,xBar0,alpha0);
end

[~,Z_mat1,nodes_col1,nodes_row1]=Z1_bound_new(A1,xBar1,alpha1);
new_Z1.mat=Z_mat1;
new_Z1.nodes_col=nodes_col1;
new_Z1.nodes_row=nodes_row1;

if nodes_col1~=nodes_col0 || nodes_row1~=nodes_row0
    [~,Z_mat0,nodes_col0,nodes_row0]=Z1_bound_new(A0,xBar0,alpha0);
    if nodes_col1~=nodes_col0 || nodes_row1~=nodes_row0
        error('not matching dimensions');
    end
end

deltax=norm_Xi_vector(xBar0-xBar1,nu)/2;
xBarhalf = (xBar1 + xBar0) ./ 2;

DDDP = Function_third_derivative(alpha0,xBarhalf,deltax);
DDP=Function_second_derivative(alpha0,xBarhalf,deltax);

% addition for the new continuation equation and linear in s scalar eqs
for i = 1:alpha1.scalar_equations.number_equations_lin
    coef_1 = extract_lin_coef(alpha1.scalar_equations, i);
    coef_0 =  extract_lin_coef(alpha0.scalar_equations, i);
    coef_Delta = coef_1-coef_0;
    DDP(i) = max(norm(vec2Xi_vec(coef_Delta,xBar0)));
end

%
% term : A_D * partial_s H_s(x_s) x_D
% 
partial_sH_sx_D = zeros(length(DDP),1);
vec_Delta = Xi_vec2vec(xBar0-xBar1);
if isintval(vec_Delta)
    partial_sH_sx_D = intval(partial_sH_sx_D);
end
for i = 1:alpha1.scalar_equations.number_equations_lin
    coef_1 = extract_lin_coef(alpha1.scalar_equations, i);
    coef_0 =  extract_lin_coef(alpha0.scalar_equations, i);
    const1 = alpha1.scalar_equations.linear_coef{3}(i);
    const0 =  alpha0.scalar_equations.linear_coef{3}(i);
    coef_Delta = coef_1-coef_0;
    const_Delta = const1-const0;
    partial_sH_sx_D(i) = abs(dot(coef_Delta,vec_Delta)+const_Delta);
end
if isintval(vec_Delta)
    partial_sH_sx_D = sup(partial_sH_sx_D);
end

% creating interval matrix for As
if isintval(A0) 
    As =  infsup (   min(inf(real(A0)),inf(real(A1))), max(sup(real(A0)),sup(real(A1))) )+...
        1i*infsup (   min(inf(imag(A0)),inf(imag(A1))), max(sup(imag(A0)),sup(imag(A1))) );
else
    As = infsup( min(real(A0),real(A1)),max(real(A0),real(A1)) )+...
        1i*infsup( min(imag(A0),imag(A1)),max(imag(A0),imag(A1)) );
end
As_ij= component_matrix_norm(As,xBar0);
Adiff_ij= component_matrix_norm(A1-A0,xBar0); % A_\Delta

sum_of_derivatives=As_ij*DDDP*deltax*deltax+2*Adiff_ij*DDP*deltax + 2*Adiff_ij*partial_sH_sx_D;

%if use_intlab
    OnesNM=(ones(N+M,1));
%else
%    OnesNM=ones(N+M,1);
%end
%addition_to_max=1/8*(2*component_matrix_norm((A1-A0)*Adagger_delta,xBar1)*OnesNM+sum_of_derivatives);
Z1DDelta=sum_of_derivatives/8;
Z1Delta=1/8*(2*component_matrix_norm((A1-A0)*Adagger_delta,xBar1)*OnesNM);
addition_to_max=Z1Delta+Z1DDelta;

if use_intlab
    Z1=test_norm_mat(max(Z_mat0.sup,Z_mat1.sup),M,N,nodes_col0,nodes_row0);
    Z1vector=Z1+addition_to_max;
else
    Z1=test_norm_mat(max(Z_mat0,Z_mat1),M,N,nodes_col0,nodes_row0);
    Z1vector=Z1+addition_to_max;
end
Z1_s= [0*Z1,Z1DDelta+Z1Delta,0*Z1, Z1];

% error handeling 
if any(Z1vector>1)
    fprintf('Z1_cont computed,\n');fprintf(' %d\n',Z1vector);
    error('Z1_cont is bigger than 1, no interval found')
elseif talkative>1
    fprintf('Z1_cont computed, %d\n',Z1vector);
elseif talkative>0
    fprintf('Z1_cont computed, %d\n',max(Z1vector));
end

