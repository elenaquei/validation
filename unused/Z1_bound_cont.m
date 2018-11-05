function [Z1vector,new_Z1,Z1_s]=Z1_bound_cont(A0,A1,xBar0,xBar1,alpha,coefs_linear0,coefs_linear1,Adagger_delta,old_Z1)
% function [Z1vector,new_Z1,Z1_s]=Z1_bound_cont(A0,A1,xBar0,xBar1,alpha,coefs_linear0,coefs_linear1,Adagger_delta,old_Z1)
% 
% In this second function to calculate the Z1 bound, another mathematical
% approach is applied. In particular, we refer to the new version of the
% attached pdf for further comments on this method. The overall idea s that
% is not interesting to split the product A[Adagger-DH]c into the norm of
% A multiplied by the norm of the rest, therefore here we will apply a
% method based on more calculations but also bringing to a
% better bound. Let us note that the previous approach was working in the
% case of van der Pol solutions but not for Lorentz. A possibility would be
% to, in each case, take the best bound between the two (this is in particular
% possible if not using intval, since the computations are relatively fast).
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
global Display

N=xBar0.size_vector;
M=xBar0.size_scalar;

% A0(1,:)=0; A1(1,:)=0;

if nargin==9 && ~isempty(old_Z1)
    Z_mat0=old_Z1.mat;
    nodes_col0=old_Z1.nodes_col;
    nodes_row0=old_Z1.nodes_row;
else
    [~,Z_mat0,nodes_col0,nodes_row0]=Z1_bound_III(A0,xBar0,alpha,coefs_linear0);
end

[~,Z_mat1,nodes_col1,nodes_row1]=Z1_bound_III(A1,xBar1,alpha,coefs_linear1);
new_Z1.mat=Z_mat1;
new_Z1.nodes_col=nodes_col1;
new_Z1.nodes_row=nodes_row1;

if nodes_col1~=nodes_col0 || nodes_row1~=nodes_row0
    [~,Z_mat0,nodes_col0,nodes_row0]=Z1_bound_III(A0,xBar0,alpha,coefs_linear0);
    if nodes_col1~=nodes_col0 || nodes_row1~=nodes_row0
        error('not matching dimensions');
    end
end

deltax=norm_Xi_vector(xBar0-xBar1,nu)/2;
xBarhalf = (xBar1 + xBar0) ./ 2;
DDDP = Function_third_derivative(xBarhalf,alpha,deltax);
DDP=Function_second_derivative(xBarhalf,alpha,deltax);

% creating interval matrix for As
if use_intlab 
    As =  infsup (   min(inf(real(A0)),inf(real(A1))), max(sup(real(A0)),sup(real(A1))) )+...
        1i*infsup (   min(inf(imag(A0)),inf(imag(A1))), max(sup(imag(A0)),sup(imag(A1))) );
else
    As = infsup( min(real(A0),real(A1)),max(real(A0),real(A1)) )+...
        1i*infsup( min(imag(A0),imag(A1)),max(imag(A0),imag(A1)) );
end
As_ij= component_matrix_norm(As,xBar0);
Adiff_ij= component_matrix_norm(A1-A0,xBar0); % A_\Delta

sum_of_derivatives=As_ij*DDDP*deltax*deltax+2*Adiff_ij*DDP*deltax;

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
elseif Display
    fprintf('Z1_cont computed, %d\n',Z1vector);
end

