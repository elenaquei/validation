function D_big = reshape_D(D,size_scalar, size_vector,old_nodes,new_nodes)
% function D_big = reshape_D(D,size_scalar, size_vector,old_nodes,new_nodes)
%
% return a bigger D matrix than the starting one.
% This is suppose to be just a normal numerical matrix that we want to
% extend. The additional elements will just be zeros ( no diagonal items)
%
%
global use_intlab

M=size_scalar;
N=size_vector;
F_D_new=new_nodes*2+1;
F=old_nodes*2+1;
diff=new_nodes-old_nodes;
%alpha1=xBar_scalar1;



if diff<0
    diff=abs(diff);
%    error('This function can just return A_out bigger than A_in');
    A_small=zeros(M+N*F_D_new);
    if use_intlab
    A_small=intval(A_small);
    end
    A_small(1:M,1:M)=D(1:M,1:M); % scalar-scalar part
    for i=1:N
        for j=1:N
            A_small(M+1+(i-1)*F_D_new:M+i*F_D_new,M+1+(j-1)*F_D_new:M+j*F_D_new)=D(M+1+(i-1)*F+diff:M+i*F-diff,M+1+(j-1)*F+diff:M+j*F-diff);
        end
        A_small(M+1+(i-1)*F_D_new:M+i*F_D_new,1:M)=D(M+1+(i-1)*F+diff:M+i*F-diff,1:M);
        A_small(1:M,M+1+(i-1)*F_D_new:M+i*F_D_new)=D(1:M,M+1+(i-1)*F+diff:M+i*F-diff);
    end
    D_big=A_small;
    return
end
if diff==0
    D_big=D;
    return
end
D_big=zeros(M+N*F_D_new);
if use_intlab
    D_big=intval(D_big);
end
D_big(1:M,1:M)=D(1:M,1:M); % scalar-scalar part

for i=0:N-1
    if ~use_intlab
    D_big(M+i*F_D_new+(1:F_D_new),1:M)=...
        [zeros(diff,M);D(M+i*F+(1:F),1:M);zeros(diff,M)];
    D_big(1:M,M+i*F_D_new+(1:F_D_new))=...
        [zeros(M,diff),D(1:M,M+i*F+(1:F)),zeros(M,diff)];
    else
        
    D_big(M+i*F_D_new+(1:F_D_new),1:M)=...
        [intval(zeros(diff,M));D(M+i*F+(1:F),1:M);intval(zeros(diff,M))];
    D_big(1:M,M+i*F_D_new+(1:F_D_new))=...
        [intval(zeros(M,diff)),D(1:M,M+i*F+(1:F)),intval(zeros(M,diff))];
    end
    for j=0:N-1
        A_ij=D(M+i*F+(1:F),M+j*F+(1:F));
        Zero1=zeros(diff);
        Zero2=zeros(F,diff);
        if use_intlab
            Zero1=intval(Zero1);
            Zero2=intval(Zero2);
        end
        if i==j
            if ~use_intlab
                K=(old_nodes+1):new_nodes;
                Diag1=diag(1i./(flip(K)));
                Diag2=diag(1i./(-K));
            else
                K=intval((old_nodes+1):new_nodes);
                Diag1=diag(intval(1i)./(flip(K)));
                Diag2=diag(intval(1i)./(-K));
            end
            D_big(M+i*F_D_new+(1:F_D_new),...
                M+j*F_D_new+(1:F_D_new))=...
                [0*Diag1,  Zero2',  Zero1;
                Zero2,   A_ij,    Zero2;
                Zero1,   Zero2',  0*Diag2];
        else
            %            if ~use_intlab
            D_big(M+i*F_D_new+(1:F_D_new),...
                M+j*F_D_new+(1:F_D_new))=...
                [Zero1,  Zero2',  Zero1;
                Zero2,   A_ij,    Zero2;
                Zero1,   Zero2',  Zero1];
            
        end 
    end
end
