function [Z1vector,Z_mat,nodes_col,nodes_row]=Z1_bound_III(A,xBar,alpha,coefs_linear)
% function Z1vector=Z1_bound_II(A,xBar,alpha,coefs_linear)
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


global use_intlab
global Display

N=xBar.size_vector; % N dimension of vector field 
M=xBar.size_scalar; % M number of scalar equations
F=2*xBar.nodes+1;   % F= (2n+1) number of nodes 
Deg=alpha.deg_vector; % Deg degree of the vector field
nodes=xBar.nodes;   
% tic;
DFm=Function_derivative_tilde(xBar,alpha,coefs_linear,1); % derivative of the vector field (without ik diagonal)
% T1=toc;
% tic
% DFm2=Function_derivative_tilde2(xBar,alpha,coefs_linear,1);
% T2=toc;
% if T1<T2
%     fprintf('tilde1 faster')
% else
%     fprintf('tilde2 faster')
%     disp(T2-T1);
%     disp(T1)
% end
%if DFm2~=DFm
%    error('major?')
%end
%DFm=mat2Xi_mat(DHm,M,N,nodes*(Deg-1));

%reshaping of the derivative, adding some nodes with zero value
if ~use_intlab
    DFm.scalar_vector=cat(3,zeros(M,N,nodes*(Deg-1)+1),DFm.scalar_vector,zeros(M,N,nodes*(Deg-1)+1));
    DFm.vector_scalar=cat(2,zeros(N,nodes*(Deg-1)+1,M),DFm.vector_scalar,zeros(N,nodes*(Deg-1)+1,M));
    DFm.vector_vector=cat(3,zeros(N,N,nodes*(Deg-1)+1),DFm.vector_vector,zeros(N,N,nodes*(Deg-1)+1));
else
    DFm.scalar_vector=intvalCAT(3,intvalCAT(3,zeros(M,N,nodes*(Deg-1)+1),DFm.scalar_vector),zeros(M,N,nodes*(Deg-1)+1));
    DFm.vector_scalar=intvalCAT(2,intvalCAT(2,zeros(N,nodes*(Deg-1)+1,M),DFm.vector_scalar),zeros(N,nodes*(Deg-1)+1,M));
    DFm.vector_vector=intvalCAT(3,intvalCAT(3,zeros(N,N,nodes*(Deg-1)+1),DFm.vector_vector),zeros(N,N,nodes*(Deg-1)+1));
end

DFm.nodes=(size(DFm.vector_scalar,2)-1)/2;%nodes+2*(nodes*(Deg-1))+1;%??

DHm=Xi_mat2mat(DFm);

nodes_new=DFm.nodes;%nodes+2*(nodes*(Deg-1))+1;%??2*(Deg-1)*nodes+1;
F_D_new=2*nodes_new+1;
diff=nodes_new-nodes;

screen=ones(size(DHm)); % this function will screen out the unwanted elements from DHm
screen(1:M,:)=0;

for i=0:N-1
    screen(M+i*F_D_new+diff+(1:F),...
        1:M)=0;
end

for i=0:N-1
    for j=0:N-1
        screen(M+i*F_D_new+diff+(1:F),...
            M+j*F_D_new+diff+(1:F))=0;
    end
end
if use_intlab
    screen=intval(screen);
end

DHm=screen.*DHm;

% at this point, the difference in size between DHm and A should be of 
% N*(Deg-2)*nodes. in order to multiply them correctly, we have to change
% the dimensions of A (remembering that on the diagonal we have the
% elements 1/(k alpha1 2pi)

% debugging check
if size(DHm,1)~=size(DHm,2)
    error('something got wrong :S');
elseif size(DHm,1) ~=  M+ N*F_D_new
    warning('probably something got wrong :-/');
end

A_big=extend_approximate_inverse(A,M,N, xBar.nodes,nodes_new);

% of DHm just the central columns need to be taken, getting to a matrix of
% the size 2(Deg-1)2nodes+3 X (Deg-1)2nodes+3
DHm_new=zeros(size(DHm,1), M+N*((Deg-1)*2*nodes+3));

if use_intlab
    DHm_new=intval(DHm_new);
end

DHm_new(:,1:M)=DHm(:,1:M);
col_new=(Deg-1)*2*nodes+3;
col_old=(size(DHm,2)-M)/N;
diff=(col_old-col_new)/2;

for col=1:N
    DHm_new(:,M+(col-1)*col_new+(1:col_new))=DHm(:,M+(col-1)*col_old+(diff+1:col_old-diff));
end

Z_mat=A_big*DHm_new;

nodes_col=(Deg-1)*nodes+1;
nodes_row=((size(Z_mat,1)-M)/N-1)/2;
Z1vector=test_norm_mat(Z_mat,M,N,nodes_col,nodes_row);

% error handeling 
if any(Z1vector>1) && Display
%    warning('Z1>1');
    fprintf('Z1 computed, %d\n',Z1vector);
    error('Z1 is bigger than 1, no interval found')
elseif Display
   fprintf('    Z1 computed, %d\n',Z1vector);
end

