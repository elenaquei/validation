classdef Xi_matrix %still to work on!!
    properties
        size_scalar
        size_vector
        scalar_scalar
        scalar_vector
        vector_scalar
        vector_vector
        nodes
    end
    methods
        % CONSTRUCTOR
        function MAT=Xi_matrix(A, scal_vec,vec_scal,vec_vec)
            % function MAT=Xi_matrix(A, scal_vec,vec_scal,vec_vec)
            % INPUT
            % A         if just one input, A is a Xi_vector (or Xi_matrix),
            % else it is a matrix, the scalar-scalar matrix, 2-d matrix
            % scal_vec  is the scalar-vector part of the Xi-matrix  3-d
            % matrix, size M,N,F
            % vec_scal  is the vector-scalar part of the Xi-matrix 3-d
            % matrix, size N,F,M
            % vec_vec   is a three-dimensional tensor for the vector-vector
            % part of the Xi-matrix of size [N,N,F]
            % NB: By the notation scalar_vector, it is meant the component of
            % the matrix that deals with vectors and returns scalars. In
            % the same way, the component vector_scalar takes as input a
            % scalar and returns a vector.
            if nargin==1
                % A is Xi-vector or Xi_matrix
                MAT.size_scalar=A.size_scalar;
                MAT.size_vector=A.size_vector;
                MAT.nodes=A.nodes;
                MAT.scalar_scalar=zeros(MAT.size_scalar);
                MAT.scalar_vector=zeros(MAT.size_scalar,MAT.size_vector,2*MAT.nodes+1);
                MAT.vector_scalar=zeros(MAT.size_vector,2*MAT.nodes+1,MAT.size_scalar);
                MAT.vector_vector=zeros(MAT.size_vector,MAT.size_vector,2*MAT.nodes+1);
            elseif nargin==2
                MAT.size_scalar=A;
                MAT.size_vector=scal_vec;
            elseif nargin==3
                MAT.size_scalar=A;
                MAT.size_vector=scal_vec;
                MAT.nodes=vec_scal;
                MAT.scalar_scalar=zeros(MAT.size_scalar);
                MAT.scalar_vector=zeros(MAT.size_scalar,MAT.size_vector,2*MAT.nodes+1);
                MAT.vector_scalar=zeros(MAT.size_vector,2*MAT.nodes+1,MAT.size_scalar);
                MAT.vector_vector=zeros(MAT.size_vector,MAT.size_vector,2*MAT.nodes+1);
            elseif nargin>1
                % A is a matrix
                if size(A,1)~=size(A,2) 
                    error('A must be square')
                elseif size(A,1) ~=size(scal_vec,1)  || size(A,1)~=size(vec_scal,3)
                    error('Dimensions must agree')
                elseif size(vec_scal,2)~=size(vec_vec,2) || size(vec_scal,1)~=size(vec_vec,1)
                    error('Number of nodes must be constant')
                elseif size(scal_vec,2)~=size(vec_vec,1) || size(scal_vec,3)~=size(vec_vec,2)
                    error('Number of nodes must be constant')
                elseif size(vec_vec,1)~=size(vec_vec,3) || size(vec_vec,2)~=size(vec_vec,4)
                    error('Dimensions must agree')
                end
                MAT.size_scalar=size(A,1);
                MAT.size_vector=size(vec_vec,1);
                MAT.nodes=(size(vec_vec,2)-1)/2+1;
                MAT.scalar_vector=scal_vec;
                MAT.vector_scalar=vec_scal;% the product in this component is to be considered as a sum over the second and third dimension
                MAT.scalar_scalar=scal_scal;
                MAT.vector_vector=vec_vec;
            end
        end
        
        % function Norm=cnorm_Xi_mat(Xi_m,nu)
        %     % component norm for Xi_matrix
        %     Norm=zeros(Xi_m.size_scalar+Xi_m.size_vector);
        %     K=-Xi_m.nodes:Xi_m.nodes;
        %     nuK=nu.^abs(K);
        %     for i=1:Xi_m.size_scalar
        %         Norm(i)=sum(abs(Xi_m.scalar,scalar),2);
        %     end
        %     for i=1:Xi_m.size_vector
        %         for j=1:Xi_m.size_vector
        %         Norm(Xi_m.size_scalar+i)=Norm(Xi_m.size_scalar+i)+...
        %             sum(abs(Xi_m.vector_vector(i,j,:)).*nuK);
        %         end
        %     end
        % end
        
        % NE
        function bool = ne(A,B)
            bool=1;
            if A.size_scalar~=B.size_scalar
                return
            elseif A.size_vector~=B.size_vector
                return
            elseif A.nodes~=B.nodes
                return
            elseif any(A.scalar_scalar~=B.scalar_scalar)
                return
            elseif any(A.scalar_vector~=B.scalar_vector)
                return
            elseif any(A.vector_scalar~=B.vector_scalar)
                return
            elseif any(A.vector_vector~=B.vector_vector)
                return
            end
            bool=0;
        end
        % NE
        
        % EQ
        function bool=eq(A,B)
            nebool=ne(A,B);
            if nebool
                bool=0;
            else
                bool=1;
            end
        end
        % EQ
        
        % INTVAL
        function A=intval(A)
            A.vector_vector=intval(A.vector_vector);
            A.vector_scalar=intval(A.vector_scalar);
            A.scalar_scalar=intval(A.scalar_scalar);
            A.scalar_vector=intval(A.scalar_vector);
        end
        % INTVAL
        
        
        % PRODUCT matrix-matrix
        function Xi_m = prod_Xi_matrix(Xi_m1,Xi_m2,flag) % NO NO NO
            % function Xi_m = prod_Xi_matrix(Xi_m1,Xi_m2,flag)
            % Xi_m1,Xi_m2     Xi_matrices with same dimensions all through
            % flag            if flag==1 (DEFAULT), the result matrix has the same
            % size as the input matrices, else it has the biggest size
            % possible, considering the tails as zeros. Eventually, the
            % only size that changes is the dimension related to the number
            % of nodes in the scal_vec or vec_vec matrices. 
            if nargin==2
                flag=1;
            end
            if Xi_m1.size_vector~=Xi_m2.size_vector 
                error('Vector dimensions must agree')
            elseif Xi_m1.size_scalar~=Xi_m2.size_scalar
                error('Scalar dimensions must agree')
            elseif Xi_m1.nodes~=Xi_m2.nodes
                if flag~=1
                    warning('Number of nodes do not agree')
                else 
                    error('Number of nodes do not agree')
                end
            end
            Xi_m=Xi_vector(Xi_m1);
            
            M=Xi_m1.size_scalar;
            N=Xi_m1.size_vector;
            F=Xi_m1.nodes;
            
            %            Xi_m.scalar_scalar [M,M]
            Xi_m.scalar_scalar=Xi_m1.scalar_scalar*Xi_m2.scalar_scalar;
            %C=zeros(M);
            for i=1:M
                for j=1:M
                    
                end
            end
            A=reshape(Xi_m1.scalar_vector,[Xi_m1.size_scalar,Xi_m1.size_vector*Xi_m1.nodes]);
            B=reshape(Xi_m2.vector_scalar,[Xi_m2.size_vector*Xi_m2.nodes,Xi_m2.size_scalar]);
            C=A*B;
            Xi_m.scalar_scalar=Xi_m.scalar_scalar+C;
            
            
            
            %            Xi_m.scalar_vector [M,N,F]
            A=Xi_m1.scalar_scalar;
            B=reshape(Xi_m2.scalar_vector,[Xi_m2.size_scalar,Xi_m2.size_vector*Xi_m2.nodes]);
            C=A*B;
            C=reshape(C,Xi_m2.size_scalar,Xi_m2.size_vector,Xi_m2.nodes);
            Xi_m.scalar_vector=C;
            % C=zeros(size_scalar,size_vector,nodes);
            A=reshape(Xi.m1.scalar_vector,[M,N*F]);
            B=reshape(Xi_m2.vector_vector,[N*F,N*F]);
            C=A*B;
            C=reshape(C,[M,N,F]);
            Xi_m.scalar_vector=Xi_m.scalar_vector+C;
            
            
            %           Xi_m.vector_scalar [N,F,M]
            A=reshape(Xi_m1.vector_scalar,[N*F,M]);
            B=Xi_m2.scalar_scalar;
            C=A*B;
            Xi_m.vector_scalar=C;
            A=reshape(Xi_m1.vector_vector,[N*F,N*F]);
            B=reshape(Xi_m2.vector_scalar,[N*F,M]);
            C=A*B;
            Xi_m.vector_scalar=Xi_m.vector_scalar+C;
            
            %            Xi_m.vector_vector
            
        end % PRODUCT matrix-matrix
        
    end
end


