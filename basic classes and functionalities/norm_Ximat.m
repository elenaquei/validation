function norm_mat=norm_Ximat(B,xBar)
% function norm_mat=norm_Ximat(B,xBar)

norm_mat=test_norm_mat(B,xBar.size_scalar,xBar.size_vector,xBar.nodes);
return

% 
% global nu
% 
% B_scal_scal=B(1:xBar.size_scalar,1:xBar.size_scalar);
% B_scal_vec=B(1:xBar.size_scalar,xBar.size_scalar+1:end);
% B_vec_scal=B(xBar.size_scalar+1:end,1:xBar.size_scalar);
% B_vec_vec=B(xBar.size_scalar+1:end,xBar.size_scalar+1:end);
% 
% norm_mat=zeros(xBar.size_scalar+xBar.size_vector,1);
% 
% norm_mat(1:xBar.size_scalar)=sum(abs(B_scal_scal),2);
% 
% for i=1:xBar.size_vector 
%     for j=1:xBar.size_scalar
%         norm_mat(j)=norm_mat(j)+...
%             InfNuNorm(B_scal_vec(j,(i-1)*(xBar.nodes*2+1)+1:i*(xBar.nodes*2+1)),nu);
%     end
% end
% 
% for j=1:xBar.size_vector
%     for i=1:xBar.size_scalar
%         norm_mat(xBar.size_scalar+j)=norm_mat(xBar.size_scalar+j)+...
%             NuNorm(B_vec_scal((j-1)*(xBar.nodes*2+1)+1:j*(xBar.nodes*2+1),i),nu);
%     end
%     for i=1:xBar.size_vector
%         Bij=B_vec_vec((j-1)*(xBar.nodes*2+1)+1:j*(xBar.nodes*2+1),...
%             (i-1)*(xBar.nodes*2+1)+1:i*(xBar.nodes*2+1));
%         norm_mat(xBar.size_scalar+j)=norm_mat(xBar.size_scalar+j)+...
%             norm_bl1(Bij,xBar.nodes,nu);
%     end
% end
% 
% 
% return
% 
% function b=norm_bl1(Bij,m,nu)
% 
% K=-m:m;
% N=nu.^abs(K);
% BB=0*K;
% 
% for n=-m:m
%     BB(n+m+1)=sum(Bij(n+m+1,:).*N)/(nu^abs(n));
% end
% 
% b=max(abs(1./N .*BB));
% 
% return