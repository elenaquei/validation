function [xBar1,x_dot1,DH1,coefs_linear1] =update(xBar0,delta,x_dot0,...
    alpha_coef,coefs_linear,maxiter,min_res)
%function [x1,x_dot1,DH1,coefs_linear1] =update(xBar0,delta,x_dot0,...
%   alpha,coefs_linear,maxiter,min_res)


% compute x_iter
x_dot0_Xi=vec2Xi_vec(x_dot0,xBar0.size_scalar,xBar0.size_vector,xBar0.nodes);
x1=xBar0+ delta* x_dot0_Xi;

% enforce reality of solution
x1=x1.symmetrise;

% add one equation for archlength continuation
coefs_linear0=coefs_linear;
coefs_linear0{1}=[-x_dot0(1:2).';[0,0]];
coefs_linear0{2}(2,:,:)=coefs_linear{2}(1,:,:);
coefs_linear0{2}(1,:,:)=reshape(-x_dot0(3:end),2*xBar0.nodes+1,xBar0.size_vector).';
coefs_linear0{3}(1) = Xi_vec2vec(x1).'*x_dot0;
coefs_linear0{3}(2)=real(-coefs_linear0{1}(2,:)*x1.scalar.'-...
    sum(sum(squeeze(coefs_linear0{2}(2,:,:)).*x1.vector,1)));


% Newton
[xBar1]=Newton_Xi(x1,alpha_coef,coefs_linear0,maxiter,min_res);
% impose symmetrisation (?)
xBar1=xBar1.symmetrise;


DH1=Function_derivative(xBar1,alpha_coef,coefs_linear,0);
% flag 0 to have non-square output

x_dot1=kernel(DH1);%null(DH1);

angle = atan( imag(x_dot1(2))/real(x_dot1(2)));
x_dot1 = exp( - 1i * angle) * x_dot1;

x_dot1=Xi_vec2vec(symmetrise(vec2Xi_vec(x_dot1,xBar0.size_scalar,xBar0.size_vector,...
xBar0.nodes)));
x_dot1=x_dot1/norm(x_dot1);

if dot(x_dot1(1:x1.size_scalar),x_dot0(1:x1.size_scalar))<0
    x_dot1=-x_dot1;
end
%if dot(x_dot1(1+x1.size_scalar:end),x_dot0(1+x1.size_scalar:end))<0
%    x_dot1(1+x1.size_scalar:end)=-x_dot1(1+x1.size_scalar:end);
%end

% theta_best=0;
% for theta=0:0.0001:2*pi
%     if norm(x_dot0-exp(1i*theta)*x_dot1)<norm(x_dot0-exp(1i*theta_best)*x_dot1)...
%             && dot(exp(1i*theta)*x_dot1,x_dot0)>0
%         theta_best=theta;
%     end
% end
% x_dot1=exp(1i*theta_best)*x_dot1;



if max(abs(DH1*x_dot1)) > 10^-6
    warning('inversion not cool')
end

if norm(x_dot0-x_dot1) > 100*delta
    warning('Norm derivatives in x0 and x1 quite big')
elseif any(norm(xBar0-xBar1) < delta/100)
    warning('Distance between consecutive points quite small, %e',norm(xBar0-xBar1))
end



coefs_linear1=coefs_linear;
coefs_linear1{1}=[-x_dot1(1:2).';[0,0]];
coefs_linear1{2}(2,:,:)=coefs_linear1{2}(1,:,:);
coefs_linear1{2}(1,:,:)=reshape(-x_dot1(3:end),2*xBar0.nodes+1,xBar0.size_vector).';
coefs_linear1{3}(1) = Xi_vec2vec(xBar1).'*x_dot1;
coefs_linear1{3}(2)=real(-coefs_linear1{1}(2,:)*xBar1.scalar.'-...
    sum(sum(squeeze(coefs_linear1{2}(2,:,:)).*xBar1.vector,1)));



DH1=[-x_dot1.';DH1];


return

end

% function xdot=symmetric(xdot,xdot0,xBar)
% 
% K=xBar.nodes;
% M=xBar.size_scalar;
% N=xBar.size_vector;
% 
% 
% 
% theta_best=0;
% 
% for theta=0:0.0001:2*pi
%     if norm(symm(exp(1i*theta)*xdot,K,M,N))<norm(symm(exp(1i*theta_best)*xdot,K,M,N))...
%             && dot(symm(exp(1i*theta_best)*xdot,K,M,N),xdot0)>0
%         theta_best=theta;
%     end
% end
% 
% xdot=exp(1i*theta_best)*xdot;
% 
% end
% 
% 
% function s=symm(x,K,M,N)
% s= [ real(x(1:M));zeros(length(x)-M,1)];
% for ii=1:N
%     indeces =  M + (2*K+1)*(ii-1) + 1:(2*K+1);
%     s(indeces)= ( x(indeces) + flip(conj(x(indeces))))/2;
% end
% end